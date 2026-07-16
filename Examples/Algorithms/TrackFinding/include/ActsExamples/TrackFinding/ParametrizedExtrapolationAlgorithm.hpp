// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>

#include <TH2D.h>

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <tbb/combinable.h>
#pragma GCC diagnostic pop

namespace Acts {
  class MagneticFieldProvider;
  }  // namespace Acts

namespace ActsExamples {

/// Example algorithm that reads/writes data from/to the event store.
class ParametrizedExtrapolationAlgorithm : public ActsExamples::IAlgorithm {
 public:
  struct Config {

    /// Input measurements collection.
    std::string inputMeasurements;

    /// Input seeds.
    std::string inputSeeds;

    /// The tracking geometry that should be used.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Output collection name.
    /// std::string output;
  };

  explicit ParametrizedExtrapolationAlgorithm(
      const Config& cfg, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Read input and copy to the output
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const override;


  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};

  ReadDataHandle<SimSeedContainer> m_inputSeeds{this, "InputSeeds"};

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  TH2D *m_cht = new TH2D("cht", "Circular Hough Transform histogram", 20000, -20000, 20000, 20000, -20000, 20000);

};

}  // namespace ActsExamples
