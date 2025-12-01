// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MultiTrackFinderError.hpp"

#include <string>

namespace {

class MultiTrackFinderErrorCategory: public std::error_category {
  public:
    // Return a short descriptive name for the category
    const char* name() const noexcept final{
      return "MultiTrackFinderError";
    }

    // Return what each enum means in text.
    std::string message(int c) const final {
      using Acts::MultiTrackFinderError;

      switch (static_cast<MultiTrackFinderError>(c)) {
        case MultiTrackFinderError::UpdateFailed:
          return "Kalman update failed";
        case MultiTrackFinderError::SmoothFailed:
          return "Kalman smooth failed";
        case MultiTrackFinderError::OutputConversionFailed:
          return "Kalman output conversion failed";
        case MultiTrackFinderError::MeasurementSelectionFailed:
          return "Measurement selection failed";
        case MultiTrackFinderError::PropagationReachesMaxSteps:
          return "Propagation reaches max steps before track finding is "
                 "finished";
        default:
          return "unknown";
      }
    }
};

} //namespace

std::error_code Acts::make_error_code(Acts::MultiTrackFinderError e) {
  static MultiTrackFinderErrorCategory c;
  return {static_cast<int>(e), c};
}