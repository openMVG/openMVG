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
    float sum = accumulate(descr, descr+128, 0.0f);
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
    float sum = accumulate(descr, descr+128, 0.0f);
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
  int firstOctave = (bDezoom == true) ? -1 : 0;
  // Number of octaves.
  int numOctaves = 6;
  // Number of scales per octave.
  int numScales = 3;
  // Max ratio of Hessian eigenvalues.
  float edgeThresh = 10.0f;
  // Min contrast.
  float peakThresh = dPeakThreshold;

  int w=I.Width(), h=I.Height();
  //Convert to float
  Image<float> If( I.GetMat().cast<float>() );


  vl_constructor();

  VlSiftFilt *filt = vl_sift_new(w,h,numOctaves,numScales,firstOctave);
  if (edgeThresh >= 0)
    vl_sift_set_edge_thresh(filt, edgeThresh);
  if (peakThresh >= 0)
    vl_sift_set_peak_thresh(filt, 255*peakThresh/numScales);

  vl_sift_process_first_octave(filt, If.data());

  vl_sift_pix descr[128];
  Descriptor<type, 128> descriptor;

  while (true) {
    vl_sift_detect(filt);

    VlSiftKeypoint const *keys  = vl_sift_get_keypoints(filt);
    int nkeys = vl_sift_get_nkeypoints(filt);

    for (int i=0;i<nkeys;++i) {
      double angles [4];
      int	nangles=vl_sift_calc_keypoint_orientations(filt,angles,keys+i);

      for (int q=0 ; q < nangles ; ++q) {
        vl_sift_calc_keypoint_descriptor(filt,descr, keys+i, angles [q]);
        SIOPointFeature fp;
        fp.x() = keys[i].x;
        fp.y() = keys[i].y;
        fp.scale() = keys[i].sigma;
        fp.orientation() = static_cast<float>(angles[q]);

        siftDescToFloat(descr, descriptor, bRootSift);
        descs.push_back(descriptor);
        feats.push_back(fp);
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
