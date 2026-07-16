// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"

#include "ActsExamples/TrackFinding/ParametrizedExtrapolationAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <math.h>
#include <TH2D.h>

#include <TFile.h>



ActsExamples::ParametrizedExtrapolationAlgorithm::ParametrizedExtrapolationAlgorithm(
    const Config& cfg, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("ParametrizedExtrapolationAlgorithm", level), m_cfg(cfg) {
  // non-optional config settings must be checked on construction.
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement collection");
  }
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds collection");
  }
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);

  m_inputSeeds.initialize(m_cfg.inputSeeds);
}

ActsExamples::ProcessCode ActsExamples::ParametrizedExtrapolationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // event-store is append-only and always returns a const reference.
  ACTS_VERBOSE("Reading HelloDataCollection " << m_cfg.inputSeeds);
  const auto& seeds = m_inputSeeds(ctx);
  ACTS_VERBOSE("Read HelloDataCollection with size " << seeds.size());

  for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    // Get the bottom space point and its reference surfce
    const auto& bottomSP = seed.sp().front();
    const auto& centreSP = seed.sp()[1];
    const auto& topSP = seed.sp().back();
    
    if (bottomSP->sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in the space point");
      continue;
    }

    // ACTS_VERBOSE(" (" << bottomSP->x() << ", " << bottomSP->y() << ", " << bottomSP->z() << "),"  <<
    //               "(" << centreSP->x() << ", " << centreSP->y() << ", " << centreSP->z() << ")" <<               
    //               "(" << topSP->x() << ", " << topSP->y() << ", " << topSP->z() << ")"
    // );

    const float u_top = topSP->x()/(pow(topSP->x(),2) + pow(topSP->y(),2));
    const float u_bottom = bottomSP->x()/(pow(bottomSP->x(),2) + pow(bottomSP->y(),2));
  
    const float v_top = topSP->y()/(pow(topSP->x(),2) + pow(topSP->y(),2));
    const float v_bottom = bottomSP->y()/(pow(bottomSP->x(),2) + pow(bottomSP->y(),2));

    float A = (v_top - v_bottom)/(u_top - u_bottom);
    float B = v_bottom - A*u_bottom;
    
    const float R = 0.5*sqrt((pow(A,2) + 1)/pow(B,2));

    float acc_size = 1000.0;
    float circ_cos, circ_sin;
    for (int itheta=0; itheta < acc_size; ++itheta){
      circ_cos = R*cos((itheta/acc_size)*(2*std::numbers::pi));
      circ_sin = R*sin((itheta/acc_size)*(2*std::numbers::pi));
      m_cht->Fill(bottomSP->x() - circ_cos, bottomSP->y() - circ_sin);
      m_cht->Fill(centreSP->x() - circ_cos,centreSP->y() - circ_sin);
      m_cht->Fill(topSP->x() - circ_cos,topSP->y() - circ_sin);
    }
    
    int x_0,y_0,z;
    long MaxBin = m_cht->GetMaximumBin();

    m_cht->GetBinXYZ(MaxBin, x_0, y_0,z);
 
    ACTS_VERBOSE("The bin having the maximum value is (" << m_cht->GetXaxis()->GetBinCenter(x_0) <<"," << m_cht->GetYaxis()->GetBinCenter(y_0) <<")");
    
    
    // Fill the histogram...
    // double max_value = h->GetMaximum();
    // std::cout << "Maximum bin content: " << max_value << std::endl;


    m_cht->Reset("");
  }

// create a copy
  //HelloDataCollection copy(in);

  // transfer the copy to the event store. this always transfers ownership
  // via r-value reference/ move construction.
  //m_writeHandle(ctx, std::move(copy));

  return ActsExamples::ProcessCode::SUCCESS;
}
