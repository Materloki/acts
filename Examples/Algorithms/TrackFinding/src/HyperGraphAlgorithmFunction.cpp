// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "ActsExamples/TrackFinding/MultiTrackFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFinding/HyperGraphAlgorithm.hpp"

#include <memory>
#include <utility>

namespace Acts {
class MagneticFieldProvider;
class TrackingGeometry;
}   // namespace Acts

namespace {

using Stepper = Acts::SympyStepper;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using MTF =
    Acts::MultiTrackFinder<Propagator, ActsExamples::TrackContainer>;

struct TrackFinderFunctionImpl
    : public ActsExamples::HyperGraphAlgorithm::TrackFinderFunction {
  MTF trackFinder;

  explicit TrackFinderFunctionImpl(MTF&& f) : trackFinder(std::move(f)) {}

  ActsExamples::HyperGraphAlgorithm::TrackFinderResult operator ()(
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::HyperGraphAlgorithm::TrackFinderOptions& options,
      ActsExamples::TrackContainer& tracks,
      ActsExamples::TrackProxy rootBranch) const override {
    return trackFinder.findTracks(initialParameters, options, tracks, 
                                  rootBranch);
  }
};

} // namespace

std::shared_ptr<ActsExamples::HyperGraphAlgorithm::TrackFinderFunction>
ActsExamples::HyperGraphAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    const Acts::Logger& logger) {
  Stepper stepper(std::move(magneticField));
  Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = false;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  MTF trackFinder(std::move(propagator), logger.cloneWithSuffix("Finder"));

  // build the track finder functions. owns the track finder object.
  return std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));
}