// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PATENTED_SIFT_SIFT_H
#define OPENMVG_PATENTED_SIFT_SIFT_H

#include <iostream>
#include <numeric>

#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"

using namespace std;
using namespace openMVG;

extern "C" {
#include "nonFree/sift/vl/sift.h"
}

namespace openMVG {

// RootSift from [1]
// [1] R. Arandjelović, A. Zisserman.
// Three things everyone should know to improve object retrieval. CVPR2012.
// -> rootsift= sqrt( sift / sum(sift) );

inline void siftDescToFloat(vl_sift_pix descr[128],
  Descriptor<float,128> & descriptor, bool brootSift = false)
{
  //rootsift= sqrt( sift / sum(sift) );
  if (brootSift)  {
    const float sum = accumulate(descr, descr+128, 0.0f);
    for (int k=0;k<128;++k)
      descriptor[k] = sqrt(descr[k]/sum);
  }
  else
    for (int k=0;k<128;++k)
      descriptor[k] = descr[k];
}


inline void siftDescToFloat(vl_sift_pix descr[128],
  Descriptor<unsigned char,128> & descriptor, bool brootSift = false)
{
  //rootsift= sqrt( sift / sum(sift) );
  if (brootSift)  {
    const float sum = accumulate(descr, descr+128, 0.0f);
    for (int k=0;k<128;++k)
      descriptor[k] = static_cast<unsigned char>(512.f*sqrt(descr[k]/sum));
  }
  else
    for (int k=0;k<128;++k)
    descriptor[k] = static_cast<unsigned char>(512.f*descr[k]);
}

template<typename type>
static bool SIFTDetector(const Image<unsigned char>& I,
  std::vector<SIOPointFeature>& feats,
  std::vector<Descriptor<type,128> >& descs,
  bool bDezoom = false,
  bool bRootSift = false,
  float dPeakThreshold = 0.04f)
{
  // First Octave Index.
  const int firstOctave = (bDezoom == true) ? -1 : 0;
  // Number of octaves.
  const int numOctaves = 6;
  // Number of scales per octave.
  const int numScales = 3;
  // Max ratio of Hessian eigenvalues.
  const float edgeThresh = 10.0f;
  // Min contrast.
  const float peakThresh = dPeakThreshold;

  const int w=I.Width(), h=I.Height();
  //Convert to float
  const Image<float> If( I.GetMat().cast<float>() );

  vl_constructor();

  VlSiftFilt *filt = vl_sift_new(w,h,numOctaves,numScales,firstOctave);
  if (edgeThresh >= 0)
    vl_sift_set_edge_thresh(filt, edgeThresh);
  if (peakThresh >= 0)
    vl_sift_set_peak_thresh(filt, 255*peakThresh/numScales);

  vl_sift_process_first_octave(filt, If.data());

  Descriptor<vl_sift_pix, 128> descr;
  Descriptor<type, 128> descriptor;

  while (true) {
    vl_sift_detect(filt);

    VlSiftKeypoint const *keys  = vl_sift_get_keypoints(filt);
    const int nkeys = vl_sift_get_nkeypoints(filt);

    // Update gradient before launching parallel extraction
    vl_sift_update_gradient(filt);

  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for private(descr, descriptor)
  #endif 
    for (int i=0;i<nkeys;++i) {
      double angles [4];
      const int	nangles=vl_sift_calc_keypoint_orientations(filt,angles,keys+i);

      for (int q=0 ; q < nangles ; ++q) {
        vl_sift_calc_keypoint_descriptor(filt,&descr[0], keys+i, angles [q]);
        const SIOPointFeature fp(keys[i].x, keys[i].y,
          keys[i].sigma, static_cast<float>(angles[q]));

        siftDescToFloat(&descr[0], descriptor, bRootSift);
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        {
          descs.push_back(descriptor);
          feats.push_back(fp);
        }
      }
    }
    if (vl_sift_process_next_octave(filt))
      break; // Last octave
  }
  vl_sift_delete(filt);

  vl_destructor();

  return true;
}


} // namespace openMVG

#endif // OPENMVG_PATENTED_SIFT_SIFT_H
