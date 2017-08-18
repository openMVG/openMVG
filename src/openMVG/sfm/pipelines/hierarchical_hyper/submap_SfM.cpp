// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_SfM.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/hypercluster.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <iomanip>


#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {


SubmapSfMReconstructionEngine::SubmapSfMReconstructionEngine(const HsfmSubmap & submap, const tracks::STLMAPTracks & map_tracks, const std::string & soutDirectory, const std::string & sloggingFile)
  : SequentialSfMReconstructionEngine(submap.sfm_data, soutDirectory, sloggingFile)
{
  map_tracks_ = map_tracks;
}

SubmapSfMReconstructionEngine::~SubmapSfMReconstructionEngine()
{}

// we override this method with an empty function since it is
// the only thing that we don't use from SequentialSfMReconstructionEngine::Process.
// In this class, map tracks is defined at construction
bool SubmapSfMReconstructionEngine::InitLandmarkTracks()
{
  return true;
}

double SubmapSfMReconstructionEngine::ComputeResidualsHistogramRedefined(Histogram<double> * histo)
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(sfm_data_.structure.size());
  for(Landmarks::const_iterator iterTracks = sfm_data_.GetLandmarks().begin();
      iterTracks != sfm_data_.GetLandmarks().end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data_.GetViews().find(itObs->first)->second.get();
      const geometry::Pose3 pose = sfm_data_.GetPoseOrDie(view);
      const std::shared_ptr<cameras::IntrinsicBase> intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      vec_residuals.push_back( fabs(residual(0)) );
      vec_residuals.push_back( fabs(residual(1)) );
    }
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.begin(), vec_residuals.end());
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl
      << "SequentialSfMReconstructionEngine::ComputeResidualsMSE." << "\n"
      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}

} // namespace sfm
} // namespace openMVG
