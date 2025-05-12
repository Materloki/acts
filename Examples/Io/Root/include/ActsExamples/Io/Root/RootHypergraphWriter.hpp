// This file is part of the Acts project.
//
// Copyright (C) 2019-2025 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include <TMatrixD.h>

class TFile;
class TTree;

namespace ActsExamples {


/// @class HypergraphWriter
///
/// Write out the information (including number of measurements, outliers, holes
/// etc., fitted track parameters and corresponding majority truth particle
/// info) of the reconstructed tracks into a TTree.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to all reconstructed tracks in one
/// single event. The event number is part of the written data.
///
/// A common file can be provided for the writer to attach his TTree, this is
/// done by setting the Config::rootFile pointer to an existing file.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootHypergraphWriter final: public WriterT<ConstTrackContainer>{
  public:
    struct Config {
    /// Input (fitted) tracks collection
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Input collection to map measured hits to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// Output filename.
    std::string filePath = "hypergraph.root";
    /// Name of the output tree.
    std::string treeName = "hypergraph";
    /// File access mode.
    std::string fileMode = "RECREATE";   
    };

    /// Constructor
    ///
    /// @param config Configuration struct
    /// @param level Message level declaration
    RootHypergraphWriter(const Config& Config, Acts::Logging::Level level);
    ~RootHypergraphWriter() override;

    /// End-of-run hook
    ProcessCode finalize() override;

    /// Get readonly access to the config parameters
    const Config& config() const { return m_cfg; }

  protected:
    /// @brief Write method called by the base class
    /// @param [in] ctx is the algorithm context for event information
    /// @param [in] tracks are what to be written out
    ProcessCode writeT(const AlgorithmContext& ctx,
                      const ConstTrackContainer& tracks) override;    

  private:
  /// The Config class
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<HitSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};
  /// Mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;
  /// The output file
  TFile* m_outputFile{nullptr};
  /// The output tree
  TTree* m_outputTree{nullptr};
  /// The event number
  std::uint32_t m_eventNr{0};
  /// The track number in event
  std::vector<std::uint32_t>  m_trackNr;

  /// The track number in event
  std::vector<std::uint32_t> m_hitId;  

  /// The track number in event
  std::uint32_t m_nStates;  
  std::uint32_t m_nMeasurements;  
  std::uint32_t m_nOutliers;  
  std::uint32_t m_nHoles;
  std::uint32_t m_nSharedHits;  
  std::uint32_t m_chi2Sum;
  std::uint32_t m_NDF;





};

} // namespace ActsExamples  