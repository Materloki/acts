// This file is part of the Acts project.
//
// Copyright (C) 2019-2025 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootHypergraphWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;


namespace ActsExamples {

RootHypergraphWriter::RootHypergraphWriter(
    const RootHypergraphWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputTracks, "RootHypergraphWriter", level),
      m_cfg(config){

   // tracks collection name is already checked by base ctor
   if (m_cfg.inputParticles.empty()) {
      throw std::invalid_argument("Missing particles input collection");
   }
   if (m_cfg.inputTrackParticleMatching.empty()) {
      throw std::invalid_argument("Missing input track particles matching");
   }
   if (m_cfg.filePath.empty()) {
      throw std::invalid_argument("Missing output filename");
   }
   if (m_cfg.treeName.empty()) {
      throw std::invalid_argument("Missing tree name");
   }


   m_inputParticles.initialize(m_cfg.inputParticles);
   m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
   m_inputSimHits.initialize(m_cfg.inputSimHits);
   m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

   // Setup ROOT I/O
   auto path = m_cfg.filePath;
   m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
   if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path + "'");
   }
   m_outputFile->cd();
   m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
   if (m_outputTree == nullptr) {
      throw std::bad_alloc();
   }
   // I/O parameters
   m_outputTree->Branch("event_nr", &m_eventNr);
   m_outputTree->Branch("track_nr", &m_trackNr);

   m_outputTree->Branch("hitId", &m_hitId);

   m_outputTree->Branch("nStates", &m_nStates);
   m_outputTree->Branch("nMeasurements", &m_nMeasurements);
   m_outputTree->Branch("nOutliers", &m_nOutliers);  
   m_outputTree->Branch("nHoles", &m_nHoles);
   m_outputTree->Branch("nSharedHits", &m_nSharedHits);
   m_outputTree->Branch("chi2", &m_chi2Sum);
   m_outputTree->Branch("DoF", &m_NDF);
}

RootHypergraphWriter::~RootHypergraphWriter() {
   m_outputFile->Close();
}

ProcessCode RootHypergraphWriter::finalize() {
   m_outputFile->cd();
   m_outputTree->Write();
   m_outputFile->Close();

   ACTS_INFO("Wrote parameters of tracks to tree '" << m_cfg.treeName << "' in '"
                                                     << m_cfg.filePath << "'");
  
   return ProcessCode::SUCCESS;   
}

ProcessCode RootHypergraphWriter::writeT(const AlgorithmContext& ctx,
                                         const ConstTrackContainer& tracks){

   // Read additional input collections
   const auto& particles = m_inputParticles(ctx);
   const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);  
   const auto& simHits = m_inputSimHits(ctx);
   const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);
   // For each particle within a track, how many hits did it contribute
   std::vector<ParticleHitCount> particleHitCounts;

   // Exclusive access to the tree while writing
   std::lock_guard<std::mutex> lock(m_writeMutex);

   // Get the event number
   m_eventNr = ctx.eventNumber;

   for (const auto& track : tracks) {
      for (const auto& state : track.trackStatesReversed()) {
         m_nStates = track.nTrackStates();
         m_nMeasurements = track.nMeasurements();
         m_nOutliers = track.nOutliers();
         m_nHoles = track.nHoles();
         m_nSharedHits = track.nSharedHits();
         m_chi2Sum = track.chi2();
         m_NDF = track.nDoF();
         if(state.hasUncalibratedSourceLink()){
            auto sl =
               state.getUncalibratedSourceLink().template get<IndexSourceLink>();
            const auto hitIdx = sl.index();
            auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
            if(!indices.empty()){
               m_hitId.push_back(indices.begin()->second);
               m_trackNr.push_back(track.index());
            } 
         }

      }
      // fill the variables
      m_outputTree->Fill();   
      m_trackNr.clear();
      m_hitId.clear();
   }



   return ProcessCode::SUCCESS;
}
} // namespace ActsExamples