// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GaussianTrackDensity3D.hpp"

#include "Acts/Utilities/Logger.hpp"

#include "Acts/Vertexing/VertexingError.hpp"
#include <math.h>
#include <stdio.h>

namespace Acts {

Result<std::optional<std::pair<Vector3, double>>>
Acts::GaussianTrackDensity3D::globalMaximumWithWidth(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto result = addTracks(state, trackList);
  if (!result.ok()) {
    return result.error();
  }

  Vector3 maxPosition = {0., 0., 0.};
  double maxDensity = 0.;
  SquareMatrix3 maxSecondDerivative;


  for (const auto& track : state.trackEntries) {
    Vector3 trial_pos = {track.pos(0), track.pos(1), track.pos(2)};

    auto [density, firstDerivative, secondDerivative] =
        trackDensityAndDerivatives(state, trial_pos);
    if ( density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trial_pos, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);

    trial_pos += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) =
        trackDensityAndDerivatives(state, trial_pos);

    if ( density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trial_pos, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
    trial_pos += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) =
        trackDensityAndDerivatives(state, trial_pos);
    if ( density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trial_pos, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
  }

  if (maxSecondDerivative.determinant() == 0.) {
    return std::nullopt;
  }

  //return std::pair{maxPosition, std::sqrt(-(maxDensity / maxSecondDerivative))};
  return std::pair{maxPosition, 0};
}

Result<void> Acts::GaussianTrackDensity3D::globalMaximum(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto maxRes = globalMaximumWithWidth(state, trackList);

  if (!maxRes.ok()) {
    return maxRes.error();
  }
  
  const auto& maxOpt = *maxRes;
  /*
  if (!maxOpt.has_value()) {
    return std::nullopt;
  }
   */
  return Result<void>::success();
}

Result<void> Acts::GaussianTrackDensity3D::addTracks(
    State& state, const std::vector<InputTrack>& trackList) const {
  for (auto trk : trackList) {
    const BoundTrackParameters& boundParams = m_cfg.extractParameters(trk);


    // Get required track parameters
    Vector3 pos;

    pos(0) = boundParams.parameters()[BoundIndices::eBoundLoc0];
    pos(1) = boundParams.parameters()[BoundIndices::eBoundLoc1];
    pos(2) = boundParams.parameters()[BoundIndices::eBoundTime];

    // Get track covariance
    if (!boundParams.covariance().has_value()) {
      return VertexingError::NoCovariance;
    }
    const auto perigeeCov = *(boundParams.covariance());
    const double covDD =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0);
    const double covZZ =
        perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1);
    const double covTT =
            perigeeCov(BoundIndices::eBoundTime, BoundIndices::eBoundTime);
    const double covDZ =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1);
    const double covDT =
            perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundTime);
    const double covZT =
            perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundTime);
    
    SquareMatrix3 covMatrix;
    covMatrix(0,0) = covDD;
    covMatrix(1,1) = covZZ;
    covMatrix(2,2) = covTT;
    covMatrix(0,1) = covDZ; covMatrix(1,0) = covDZ;
    covMatrix(0,2) = covDT; covMatrix(2,0) = covDT;
    covMatrix(1,2) = covZT; covMatrix(2,1) = covZT;

    state.trackEntries.emplace_back(pos, covMatrix.inverse(), covMatrix.determinant(),
                                    0, 0);
  }
  return Result<void>::success();
}

std::tuple<double, Vector3, SquareMatrix3>
Acts::GaussianTrackDensity3D::trackDensityAndDerivatives(State& state,
                                                       Vector3 pos) const {
  GaussianTrackDensityStore densityResult(pos);
  for (const auto& trackEntry : state.trackEntries) {
    densityResult.addTrackToDensity(trackEntry);
  }
  return densityResult.densityAndDerivatives();
}

std::tuple<Vector3, double, SquareMatrix3> Acts::GaussianTrackDensity3D::updateMaximum(
    Vector3 newPos, double newValue, SquareMatrix3 newHermitiam, Vector3 maxPos,
    double maxValue, SquareMatrix3 maxHermitiam) const {
  if (newValue > maxValue) {
    maxPos = newPos;
    maxValue = newValue;
    maxHermitiam = newHermitiam;
  }
  return {maxPos, maxValue, maxHermitiam};
}

Vector3 Acts::GaussianTrackDensity3D::stepSize( double density,
                                            Vector3 gradient,
                                            SquareMatrix3 hermitiam) const {
  //return (m_cfg.isGaussianShaped ? (y * dy) / (dy * dy - y * ddy) : -dy / ddy);
  if(m_cfg.isGaussianShaped){

    SquareMatrix3 denominator = gradient*gradient.transpose() - hermitiam*density;

    return(density*gradient.transpose()*denominator.inverse());

  }
  else{
    return (-gradient.transpose()*hermitiam.inverse());
  }                                           
}

void Acts::GaussianTrackDensity3D::GaussianTrackDensityStore::addTrackToDensity(
    const TrackEntry& entry) {
  
  
  Vector3 delta_pos = m_pos -entry.pos; 

  double deltaW = std::exp(-0.5*delta_pos.transpose()*entry.covInv*delta_pos);
  m_density += deltaW; 

  double inv_cov_a = entry.covInv(0,0);
  double inv_cov_db = entry.covInv(1,0) + entry.covInv(0,1);
  double inv_cov_e = entry.covInv(1,1);
  double inv_cov_fh = entry.covInv(1,2) + entry.covInv(2,1);
  double inv_cov_i = entry.covInv(2,2);
  double inv_cov_gc = entry.covInv(0,2) + entry.covInv(2,0);


  double r_partial = -entry.pos(0)*2*inv_cov_a + entry.pos(1)*inv_cov_db + entry.pos(2)*inv_cov_gc;
  double z_partial = -entry.pos(0)*inv_cov_db + entry.pos(1)*2*inv_cov_e + entry.pos(2)*inv_cov_fh;
  double t_partial = -entry.pos(0)*inv_cov_gc + entry.pos(1)*2*inv_cov_fh + entry.pos(2)*2*inv_cov_i;
  
  m_grad(0) += -0.5*deltaW*z_partial;
  m_grad(1) += -0.5*deltaW*t_partial;    
  m_grad(2) += -0.5*deltaW*r_partial;


  m_hermitian(0,1) +=  - 0.25*deltaW*(-2*inv_cov_db + z_partial);
  m_hermitian(0,2) +=  - 0.25*deltaW*(-2*inv_cov_gc + t_partial);

  m_hermitian(1,0) +=  - 0.25*deltaW*(-2*inv_cov_db + r_partial);
  m_hermitian(1,2) +=  - 0.25*deltaW*(-2*inv_cov_fh + t_partial);
  
  m_hermitian(2,0) +=  - 0.25*deltaW*(-2*inv_cov_gc + r_partial);
  m_hermitian(2,1) +=  - 0.25*deltaW*(-2*inv_cov_fh + z_partial);

  m_hermitian(0,0) +=  - 0.25*deltaW*(-4*inv_cov_a + r_partial);
  m_hermitian(1,1) +=  - 0.25*deltaW*(-4*inv_cov_e + z_partial);
  m_hermitian(2,2) +=  - 0.25*deltaW*(-4*inv_cov_i + t_partial);
  

}

}  // namespace Acts
