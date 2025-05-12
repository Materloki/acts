// This file is part of the Acts project.
//
// Copyright (C) 2020-2025 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/HyperGraphAlgorithm.hpp"
#include "ActsExamples/TrackFinding/MultiTrackNavigator.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <boost/functional/hash.hpp>

// Specialize std::hash for SeedIdentifier
// This is required to use SeedIdentifier as a key in an `std::unordered_map`.
template <class T, std::size_t N>
struct std::hash<std::array<T, N>> {
  std::size_t operator()(const std::array<T, N>& array) const {
    std::hash<T> hasher;
    std::size_t result = 0;
    for (auto&& element : array) {
      boost::hash_combine(result, hasher(element));
    }
    return result;
  }
};

namespace ActsExamples {

namespace {

class MeasurementSelector {
 public:
  using Traj = Acts::VectorMultiTrajectory;

  explicit MeasurementSelector(Acts::MeasurementSelector selector)
      : m_selector(std::move(selector)) {}

  void setSeed(const std::optional<SimSeed>& seed) { m_seed = seed; }

  Acts::Result<std::pair<std::vector<Traj::TrackStateProxy>::iterator,
                         std::vector<Traj::TrackStateProxy>::iterator>>
  select(std::vector<Traj::TrackStateProxy>& candidates, bool& isOutlier,
         const Acts::Logger& logger) const {
    if (m_seed.has_value()) {
      std::vector<Traj::TrackStateProxy> newCandidates;

      for (const auto& candidate : candidates) {
        if (isSeedCandidate(candidate)) {
          newCandidates.push_back(candidate);
        }
      }

      if (!newCandidates.empty()) {
        candidates = std::move(newCandidates);
      }
    }

    return m_selector.select<Acts::VectorMultiTrajectory>(candidates, isOutlier,
                                                          logger);
  }

 private:
  Acts::MeasurementSelector m_selector;
  std::optional<SimSeed> m_seed;

  bool isSeedCandidate(const Traj::TrackStateProxy& candidate) const {
    assert(candidate.hasUncalibratedSourceLink());

    const Acts::SourceLink& sourceLink = candidate.getUncalibratedSourceLink();

    for (const auto& sp : m_seed->sp()) {
      for (const auto& sl : sp->sourceLinks()) {
        if (sourceLink.get<IndexSourceLink>() == sl.get<IndexSourceLink>()) {
          return true;
        }
      }
    }

    return false;
  }
};

/// Source link indices of the bottom, middle, top measurements.
/// In case of strip seeds only the first source link of the pair is used.
using SeedIdentifier = std::array<Index, 3>;

/// Build a seed identifier from a seed.
///
/// @param seed The seed to build the identifier from.
/// @return The seed identifier.
SeedIdentifier makeSeedIdentifier(const SimSeed& seed) {
  SeedIdentifier result;

  for (const auto& [i, sp] : Acts::enumerate(seed.sp())) {
    const Acts::SourceLink& firstSourceLink = sp->sourceLinks().front();
    result.at(i) = firstSourceLink.get<IndexSourceLink>().index();
  }

  return result;
}

/// Visit all possible seed identifiers of a track.
///
/// @param track The track to visit the seed identifiers of.
/// @param visitor The visitor to call for each seed identifier.
template <typename Visitor>
void visitSeedIdentifiers(const TrackProxy& track, Visitor visitor) {
  // first we collect the source link indices of the track states
  std::vector<Index> sourceLinkIndices;
  sourceLinkIndices.reserve(track.nMeasurements());
  for (const auto& trackState : track.trackStatesReversed()) {
    if (!trackState.hasUncalibratedSourceLink()) {
      continue;
    }
    const Acts::SourceLink& sourceLink = trackState.getUncalibratedSourceLink();
    sourceLinkIndices.push_back(sourceLink.get<IndexSourceLink>().index());
  }

  // then we iterate over all possible triplets and form seed identifiers
  for (std::size_t i = 0; i < sourceLinkIndices.size(); ++i) {
    for (std::size_t j = i + 1; j < sourceLinkIndices.size(); ++j) {
      for (std::size_t k = j + 1; k < sourceLinkIndices.size(); ++k) {
        // Putting them into reverse order (k, j, i) to compensate for the
        // `trackStatesReversed` above.
        visitor({sourceLinkIndices.at(k), sourceLinkIndices.at(j),
                 sourceLinkIndices.at(i)});
      }
    }
  }
}

class BranchStopper {
 public:
  using BranchStopperResult =
      Acts::MultiTrackNavigatorBranchStopperResult;

  struct BrachState {
    std::size_t nPixelHoles = 0;
    std::size_t nStripHoles = 0;
  };

  static constexpr Acts::ProxyAccessor<BrachState> branchStateAccessor =
      Acts::ProxyAccessor<BrachState>(Acts::hashString("MyBranchState"));

  mutable std::atomic<std::size_t> m_nStoppedBranches{0};

  explicit BranchStopper(const HyperGraphAlgorithm::Config& config)
      : m_cfg(config) {}

  BranchStopperResult operator()(
      const TrackContainer::TrackProxy& track,
      const TrackContainer::TrackStateProxy& trackState) const {
    if (!m_cfg.trackSelectorCfg.has_value()) {
      return BranchStopperResult::Continue;
    }

    const Acts::TrackSelector::Config* singleConfig = std::visit(
        [&](const auto& config) -> const Acts::TrackSelector::Config* {
          using T = std::decay_t<decltype(config)>;
          if constexpr (std::is_same_v<T, Acts::TrackSelector::Config>) {
            return &config;
          } else if constexpr (std::is_same_v<
                                   T, Acts::TrackSelector::EtaBinnedConfig>) {
            double theta = trackState.parameters()[Acts::eBoundTheta];
            double eta = -std::log(std::tan(0.5 * theta));
            return config.hasCuts(eta) ? &config.getCuts(eta) : nullptr;
          }
        },
        *m_cfg.trackSelectorCfg);

    if (singleConfig == nullptr) {
      ++m_nStoppedBranches;
      return BranchStopperResult::StopAndDrop;
    }

    bool tooManyHolesPS = false;
    if (!(m_cfg.pixelVolumes.empty() && m_cfg.stripVolumes.empty())) {
      auto& branchState = branchStateAccessor(track);
      // count both holes and outliers as holes for pixel/strip counts
      if (trackState.typeFlags().test(Acts::TrackStateFlag::HoleFlag) ||
          trackState.typeFlags().test(Acts::TrackStateFlag::OutlierFlag)) {
        if (m_cfg.pixelVolumes.count(
                trackState.referenceSurface().geometryId().volume()) >= 1) {
          ++branchState.nPixelHoles;
        } else if (m_cfg.stripVolumes.count(
                       trackState.referenceSurface().geometryId().volume()) >=
                   1) {
          ++branchState.nStripHoles;
        }
      }
      tooManyHolesPS = branchState.nPixelHoles > m_cfg.maxPixelHoles ||
                       branchState.nStripHoles > m_cfg.maxStripHoles;
    }

    bool enoughMeasurements =
        track.nMeasurements() >= singleConfig->minMeasurements;
    bool tooManyHoles =
        track.nHoles() > singleConfig->maxHoles || tooManyHolesPS;
    bool tooManyOutliers = track.nOutliers() > singleConfig->maxOutliers;

    if (tooManyHoles || tooManyOutliers) {
      ++m_nStoppedBranches;
      return enoughMeasurements ? BranchStopperResult::StopAndKeep
                                : BranchStopperResult::StopAndDrop;
    }

    return BranchStopperResult::Continue;
  }

 private:
  const HyperGraphAlgorithm::Config& m_cfg;
};
} // namespace

HyperGraphAlgorithm::HyperGraphAlgorithm(Config config,
                                            Acts::Logging::Level level)
    : IAlgorithm("HyperGraphAlgorithm", level), m_cfg(std::move(config)) {
    if (m_cfg.inputMeasurements.empty()) {
        throw std::invalid_argument("Missing measurements input collection");
    }
    if (m_cfg.inputSourceLinks.empty()) {
        throw std::invalid_argument("Missing source links input collection");
    }
    if (m_cfg.inputInitialTrackParameters.empty()) {
        throw std::invalid_argument(
            "Missing initial track parameters input collection");
    }
    m_inputMeasurements.initialize(m_cfg.inputMeasurements);
    m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
    m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
    m_inputSeeds.maybeInitialize(m_cfg.inputSeeds);
    m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode HyperGraphAlgorithm::execute(const AlgorithmContext& ctx) const {
    const auto& measurements = m_inputMeasurements(ctx);
    const auto& sourceLinks = m_inputSourceLinks(ctx);
    const auto& initialParameters = m_inputInitialTrackParameters(ctx);
    const SimSeedContainer* seeds = nullptr;

    if (m_inputSeeds.isInitialized()) {
        seeds = &m_inputSeeds(ctx);

        if (initialParameters.size() != seeds->size()) {
        ACTS_ERROR("Number of initial parameters and seeds do not match. "
                    << initialParameters.size() << " != " << seeds->size());
        }
    }

    // Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3{0., 0., 0.});

    PassThroughCalibrator pcalibrator;
    MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);

    MeasurementSelector measSel{
    Acts::MeasurementSelector(m_cfg.measurementSelectorCfg)};

    BranchStopper branchStopper(m_cfg);

    using Extensions = Acts::MultiTrackNavigatorExtensions<TrackContainer>;
    Extensions extensions;
    extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
    extensions.measurementSelector.connect<&MeasurementSelector::select>(
        &measSel);
    extensions.branchStopper.connect<&BranchStopper::operator()>(&branchStopper);

    IndexSourceLinkAccessor slAccessor;
    slAccessor.container = &sourceLinks;
    Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
        slAccessorDelegate;
    slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

    Acts::PropagatorPlainOptions firstPropOptions(ctx.geoContext,
                                                  ctx.magFieldContext);
    firstPropOptions.maxSteps = m_cfg.maxSteps;
    firstPropOptions.direction = Acts::Direction::Forward;

    // Set the CombinatorialKalmanFilter options
    HyperGraphOptions firstOptions(ctx.geoContext, ctx.magFieldContext,
                                    ctx.calibContext, slAccessorDelegate,
                                    extensions, firstPropOptions);
    firstOptions.targetSurface = nullptr;


    using Extrapolator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
    using ExtrapolatorOptions =
        Extrapolator::template Options<Acts::ActionList<Acts::MaterialInteractor>,
                                        Acts::AbortList<Acts::EndOfWorldReached>>;

    Extrapolator extrapolator(
        Acts::SympyStepper(m_cfg.magneticField),
        Acts::Navigator({m_cfg.trackingGeometry},
                        logger().cloneWithSuffix("Navigator")),
        logger().cloneWithSuffix("Propagator"));

    ExtrapolatorOptions extrapolationOptions(ctx.geoContext, ctx.magFieldContext);

    ACTS_DEBUG("Invoked HyperGraph builder with " << seeds->size() << " seeds");

    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

    auto trackContainerTemp = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainerTemp =
        std::make_shared<Acts::VectorMultiTrajectory>();

    TrackContainer tracks(trackContainer, trackStateContainer);
    TrackContainer tracksTemp(trackContainerTemp, trackStateContainerTemp);

    // Note that not all backends support PODs as column types
    tracks.addColumn<BranchStopper::BrachState>("MyBranchState");
    tracksTemp.addColumn<BranchStopper::BrachState>("MyBranchState");

    tracks.addColumn<unsigned int>("trackGroup");
    tracksTemp.addColumn<unsigned int>("trackGroup");
    Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");
    unsigned int nSeed = 0;

    // A map indicating whether a seed has been discovered already
    std::unordered_map<SeedIdentifier, bool> discoveredSeeds;

    auto addTrack = [&](const TrackProxy& track) {
        ++m_nFoundTracks;

        // flag seeds which are covered by the track
        visitSeedIdentifiers(track, [&](const SeedIdentifier& seedIdentifier) {
        if (auto it = discoveredSeeds.find(seedIdentifier);
            it != discoveredSeeds.end()) {
            it->second = true;
        }
        });

        if (m_trackSelector.has_value() && !m_trackSelector->isValidTrack(track)) {
        return;
        }

        ++m_nSelectedTracks;

        auto destProxy = tracks.makeTrack();
        // make sure we copy track states!
        destProxy.copyFrom(track, true);
    };

    if (seeds != nullptr && m_cfg.seedDeduplication) {
        // Index the seeds for deduplication
        for (const auto& seed : *seeds) {
        SeedIdentifier seedIdentifier = makeSeedIdentifier(seed);
        discoveredSeeds.emplace(seedIdentifier, false);
        }
    }

    for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
        m_nTotalSeeds++;
        if (seeds != nullptr) {
        const SimSeed& seed = seeds->at(iSeed);

        if (m_cfg.seedDeduplication) {
            SeedIdentifier seedIdentifier = makeSeedIdentifier(seed);
            // check if the seed has been discovered already
            if (auto it = discoveredSeeds.find(seedIdentifier);
                it != discoveredSeeds.end() && it->second) {
            m_nDeduplicatedSeeds++;
            ACTS_VERBOSE("Skipping seed " << iSeed << " due to deduplication.");
            continue;
            }
        }

        if (m_cfg.stayOnSeed) {
            measSel.setSeed(seed);
        }
        }

        // Clear trackContainerTemp and trackStateContainerTemp
        tracksTemp.clear();

        const Acts::BoundTrackParameters& InitialParameters =
            initialParameters.at(iSeed);

        ACTS_VERBOSE("iSeed " << iSeed << "\n"
                   << " Initial Parameters \n" << InitialParameters);
        auto firstResult =
            (*m_cfg.findTracks)(InitialParameters, firstOptions, tracksTemp);
        nSeed++;

        if (!firstResult.ok()) {
          m_nFailedSeeds++;
          ACTS_WARNING("Track finding failed for seed " << iSeed << " with error"
                                                        << firstResult.error());
          continue;
        }        
      auto& firstTracksForSeed = firstResult.value();
      for (auto& firstTrack : firstTracksForSeed) {
        // TODO a copy of the track should not be necessary but is the safest way
        //      with the current EDM
        // TODO a lightweight copy without copying all the track state components
        //      might be a solution
        auto trackCandidate = tracksTemp.makeTrack();
        trackCandidate.copyFrom(firstTrack, true);

        // Set the seed number, this number decrease by 1 since the seed number
        // has already been updated
        seedNumber(trackCandidate) = nSeed - 1;
        auto firstExtrapolationResult =
            Acts::extrapolateTrackToReferenceSurface(
                trackCandidate, *pSurface, extrapolator, extrapolationOptions,
                m_cfg.extrapolationStrategy, logger());
        if (!firstExtrapolationResult.ok()) {
          m_nFailedExtrapolation++;
          ACTS_ERROR("Extrapolation for seed "
                      << iSeed << " and track " << firstTrack.index()
                      << " failed with error "
                      << firstExtrapolationResult.error());
          continue;
        }
        addTrack(trackCandidate);        
      }
    }

  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");

  m_nStoppedBranches += branchStopper.m_nStoppedBranches;

  m_memoryStatistics.local().hist +=
      tracks.trackStateContainer().statistics().hist;

  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTracks(ctx, std::move(constTracks));

    
    return ProcessCode::SUCCESS;
}

ProcessCode HyperGraphAlgorithm::finalize(){
  ACTS_INFO("TrackFindingAlgorithm statistics:");
  ACTS_INFO("- total seeds: " << m_nTotalSeeds);
  ACTS_INFO("- deduplicated seeds: " << m_nDeduplicatedSeeds);
  ACTS_INFO("- failed seeds: " << m_nFailedSeeds);
  ACTS_INFO("- failed smoothing: " << m_nFailedSmoothing);
  ACTS_INFO("- failed extrapolation: " << m_nFailedExtrapolation);
  ACTS_INFO("- failure ratio seeds: " << static_cast<double>(m_nFailedSeeds) /
                                             m_nTotalSeeds);
  ACTS_INFO("- found tracks: " << m_nFoundTracks);
  ACTS_INFO("- selected tracks: " << m_nSelectedTracks);
  ACTS_INFO("- stopped branches: " << m_nStoppedBranches);

  auto memoryStatistics =
      m_memoryStatistics.combine([](const auto& a, const auto& b) {
        Acts::VectorMultiTrajectory::Statistics c;
        c.hist = a.hist + b.hist;
        return c;
      });
  std::stringstream ss;
  memoryStatistics.toStream(ss);
  ACTS_DEBUG("Track State memory statistics (averaged):\n" << ss.str());  
  return ProcessCode::SUCCESS;
}

} //namespace ActsExamples