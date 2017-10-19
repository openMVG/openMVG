// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "image_describer_akaze.hpp"

#include "openMVG/features/akaze/mldb_descriptor.hpp"
#include "openMVG/features/akaze/msurf_descriptor.hpp"
#include "openMVG/features/liop/liop_descriptor.hpp"

namespace openMVG {
namespace features {

std::unique_ptr<AKAZE_Image_describer_SURF::Regions_type>
AKAZE_Image_describer_SURF::Describe_AKAZE_SURF
(
  const image::Image<unsigned char>& image,
  const image::Image<unsigned char>* mask
)
{
  auto regions = std::unique_ptr<Regions_type>(new Regions_type);

  if (image.size() == 0)
    return regions;

  params_.options_.fDesc_factor = GetfDescFactor();

  AKAZE akaze(image, params_.options_);
  akaze.Compute_AKAZEScaleSpace();
  std::vector<AKAZEKeypoint> kpts;
  kpts.reserve(5000);
  akaze.Feature_Detection(kpts);
  akaze.Do_Subpixel_Refinement(kpts);

  // Feature masking (remove keypoints if they are masked)
  kpts.erase(std::remove_if(kpts.begin(),
                            kpts.end(),
                            [&](const AKAZEKeypoint & pt)
                            {
                              if (mask) return (*mask)(pt.y, pt.x) == 0;
                              else return false;
                            }),
             kpts.end());

  regions->Features().resize(kpts.size());
  regions->Descriptors().resize(kpts.size());

  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int i = 0; i<static_cast<int>(kpts.size()); ++i) {
    AKAZEKeypoint ptAkaze = kpts[i];

    const TEvolution& cur_slice = akaze.getSlices()[ptAkaze.class_id];

    if (bOrientation_)
      akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
    else
      ptAkaze.angle = 0.0f;

    regions->Features()[i] =
      SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

    ComputeMSURFDescriptor(cur_slice.Lx, cur_slice.Ly, ptAkaze.octave,
      regions->Features()[i],
      regions->Descriptors()[i]);
  }
  return regions;
}

std::unique_ptr<AKAZE_Image_describer_LIOP::Regions_type>
AKAZE_Image_describer_LIOP::Describe_AKAZE_LIOP
(
  const image::Image<unsigned char>& image,
  const image::Image<unsigned char>* mask
)
{
  auto regions = std::unique_ptr<Regions_type>(new Regions_type);

  if (image.size() == 0)
    return regions;

  params_.options_.fDesc_factor = GetfDescFactor();

  AKAZE akaze(image, params_.options_);
  akaze.Compute_AKAZEScaleSpace();
  std::vector<AKAZEKeypoint> kpts;
  kpts.reserve(5000);
  akaze.Feature_Detection(kpts);
  akaze.Do_Subpixel_Refinement(kpts);

  // Feature masking (remove keypoints if they are masked)
  kpts.erase(std::remove_if(kpts.begin(),
                            kpts.end(),
                            [&](const AKAZEKeypoint & pt)
                            {
                              if (mask) return (*mask)(pt.y, pt.x) == 0;
                              else return false;
                            }),
             kpts.end());

  regions->Features().resize(kpts.size());
  regions->Descriptors().resize(kpts.size());

  // Init LIOP extractor
  LIOP::Liop_Descriptor_Extractor liop_extractor;

  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int i = 0; i<static_cast<int>(kpts.size()); ++i) {
    AKAZEKeypoint ptAkaze = kpts[i];

    const TEvolution& cur_slice = akaze.getSlices()[ptAkaze.class_id];

    if (bOrientation_)
      akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
    else
      ptAkaze.angle = 0.0f;

    regions->Features()[i] =
      SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

    // Compute LIOP descriptor (do not need rotation computation, since
    //  LIOP descriptor is rotation invariant).
    // Rescale for LIOP patch extraction
    const SIOPointFeature fp =
      SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size/2.0, ptAkaze.angle);

    float desc[144];
    liop_extractor.extract(image, fp, desc);
    for (int j = 0; j<144; ++j)
      regions->Descriptors()[i][j] =
        static_cast<unsigned char>(desc[j]*255.f+.5f);
  }
  return regions;
}

std::unique_ptr<AKAZE_Image_describer_MLDB::Regions_type>
AKAZE_Image_describer_MLDB::Describe_AKAZE_MLDB
(
  const image::Image<unsigned char>& image,
  const image::Image<unsigned char>* mask
)
{
  auto regions = std::unique_ptr<Regions_type>(new Regions_type);

  if (image.size() == 0)
    return regions;

  params_.options_.fDesc_factor = GetfDescFactor();

  AKAZE akaze(image, params_.options_);
  akaze.Compute_AKAZEScaleSpace();
  std::vector<AKAZEKeypoint> kpts;
  kpts.reserve(5000);
  akaze.Feature_Detection(kpts);
  akaze.Do_Subpixel_Refinement(kpts);

  // Feature masking (remove keypoints if they are masked)
  kpts.erase(std::remove_if(kpts.begin(),
                            kpts.end(),
                            [&](const AKAZEKeypoint & pt)
                            {
                              if (mask) return (*mask)(pt.y, pt.x) == 0;
                              else return false;
                            }),
             kpts.end());

  regions->Features().resize(kpts.size());
  regions->Descriptors().resize(kpts.size());

  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int i = 0; i<static_cast<int>(kpts.size()); ++i) {
    AKAZEKeypoint ptAkaze = kpts[i];

    const TEvolution& cur_slice = akaze.getSlices()[ptAkaze.class_id];

    if (bOrientation_)
      akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
    else
      ptAkaze.angle = 0.0f;

    regions->Features()[i] =
      SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

    // Compute MLDB descriptor
    Descriptor<bool, 486> desc;
    ComputeMLDBDescriptor(cur_slice.cur, cur_slice.Lx, cur_slice.Ly,
      ptAkaze.octave, regions->Features()[i], desc);
    // convert the bool vector to the binary unsigned char array
    unsigned char* ptr = reinterpret_cast<unsigned char*>(&regions->Descriptors()[i]);
    memset(ptr, 0, regions->DescriptorLength()*sizeof(unsigned char));
    // For each byte
    for (int j = 0; j<std::ceil(486./8.); ++j, ++ptr)
    {
      // set the corresponding 8bits to the good values
      for (int iBit = 0; iBit<8 && j*8+iBit<486; ++iBit)
      {
        *ptr |= desc[j*8+iBit] << iBit;
      }
    }
  }
  return regions;
}

// static
std::unique_ptr<AKAZE_Image_describer> AKAZE_Image_describer::create
(
  const AKAZE_Image_describer::Params& params,
  bool orientation
)
{
  switch (params.eAkazeDescriptor_) {
  case AKAZE_MSURF:
    return std::unique_ptr<AKAZE_Image_describer>
        (new AKAZE_Image_describer_SURF(params, orientation));
  case AKAZE_LIOP:
    return std::unique_ptr<AKAZE_Image_describer>
        (new AKAZE_Image_describer_LIOP(params, orientation));
  case AKAZE_MLDB:
    return std::unique_ptr<AKAZE_Image_describer>
        (new AKAZE_Image_describer_MLDB(params, orientation));
  default:
    return {};
  }
}

} // namespace features
} // namespace openMVG
