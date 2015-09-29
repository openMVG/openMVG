
// Copyright (C) 2013  "Robot Vision Group, NLPR, CASIA", Zhenhua Wang,
// Bin Fan and Fuchao Wu.

// Copyright (C) 2014: Adaptation to openMVG by Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//------------------
//-- Bibliography --
//------------------
//- [1] "Local Intensity Order Pattern for Feature Description"
//- Authors: Zhenhua Wang, Bin Fan and Fuchao Wu
//- Date: 2011, ICCV, IEEE International Conference on Computer Vision

#include <openMVG/features/liop/liop_descriptor.hpp>

#include <float.h>
#include <algorithm>

namespace openMVG {
namespace features{
namespace LIOP    {

const int MAX_REGION_NUM = 10;
const int MAX_PIXEL_NUM  = 1681;
const int MAX_SAMPLE_NUM = 10;
const int LIOP_NUM = 4;
const int REGION_NUM = 6;

Liop_Descriptor_Extractor::Liop_Descriptor_Extractor()
{
  GeneratePatternMap( m_LiopPatternMap, m_LiopPosWeight, LIOP_NUM);
}

struct Pixel{
  float x, y; // position
  float f_gray; // color value
  float weight;
  int l_pattern;
};

inline bool fGrayComp(const Pixel & p1, const Pixel & p2)
{
  return p1.f_gray < p2.f_gray;
}

bool BilinearInterpolation_BorderCheck(
  float& val,
  float x, float y,
  const image::Image<float> & image,
  const image::Image<unsigned char> & flagImage)
{
  val = 0.0f;
  if(!(x >= 0 && y >= 0 && x<=image.Width()-1 && y<=image.Height()-1))
    return false;

  const int x1 = (int)x;
  const int y1 = (int)y;

  int x2, y2;
  if( x1 == image.Width()-1)
    x2  = x1;
  else
    x2 = x1+1;

  if(y1 == image.Height()-1 )
    y2 = y1;
  else
    y2 = y1+1;

  const int step = image.Width();
  const float* data = image.data();

  const int flag_step = flagImage.Width();
  const unsigned char* flag_data = flagImage.data();

  if(flag_data[y1*flag_step+x1] == 0 ||
     flag_data[y1*flag_step+x2] == 0 ||
     flag_data[y2*flag_step+x1] == 0 ||
     flag_data[y2*flag_step+x2] == 0)
    return false;

  val =
    (x2 - x) * (y2 - y) * data[y1*step+x1] +
    (x - x1) * (y2 - y) * data[y1*step+x2] +
    (x2 - x) * (y - y1) * data[y2*step+x1] +
    (x - x1) * (y - y1) * data[y2*step+x2];
  return true;
}

//non descend
void SortGray(float* dst, int* idx, float* src, int len)
{
  int i, j;

  for (i=0; i<len; ++i)
  {
    dst[i] = src[i];
    idx[i] = i;
  }

  for (i=0; i<len; ++i)
  {
    float temp = dst[i];
    int tempIdx = idx[i];
    for (j=i+1; j<len; ++j)
    {
      if (dst[j]<temp)
      {
        temp = dst[j];
        dst[j] = dst[i];
        dst[i] = temp;

        tempIdx = idx[j];
        idx[j] = idx[i];
        idx[i] = tempIdx;
      }
    }
  }
}

/// To avoid linear illumination change, we normalize the descriptor.
/// To avoid non-linear illumination change, we threshold the value
/// of each descriptor element to 'illuThresh', then normalize again.
void ThreshNorm(float* des, int dim, float illuThresh)
{
  // Normalize the descriptor, and threshold
  // value of each element to 'illuThresh'.

  float norm = 0.f;
  int i;

  for (i=0; i<dim; ++i)
  {
    norm += des[i] * des[i];
  }
  norm = sqrt(norm);
  assert(norm != 0);

  if (illuThresh <  1.f)
  {
    for (i=0; i<dim; ++i)
    {
      des[i] /= norm;

      if (des[i] > illuThresh)
      {
        des[i] = illuThresh;
      }
    }

    // Normalize again.

    norm = 0.f;

    for (i=0; i<dim; ++i)
    {
      norm += des[i] * des[i];
    }

    norm = sqrt(norm);
    assert(norm != 0);
  }

  for (i=0; i<dim; ++i)
  {
    des[i] /= norm;
  }
}

void Liop_Descriptor_Extractor::CreateLIOP_GOrder(
  const image::Image<float> & outPatch,
  const image::Image<unsigned char> & flagPatch,
  const int inRadius,
  float desc[144]) const
{
  const float* out_data = outPatch.data();
  const int out_step = outPatch.Width();
  const unsigned char* flag_data = flagPatch.data();
  const int flag_step = flagPatch.Width();

  const float inRadius2 = float(inRadius*inRadius);
  const int outRadius = outPatch.Width() / 2;

  const int lsRadius = 6;
  const float theta = 2.0f*M_PI/(float)LIOP_NUM;

  Pixel pixel[MAX_PIXEL_NUM];
  int pixelCount = 0;
  int idx[MAX_SAMPLE_NUM];
  float dst[MAX_SAMPLE_NUM];
  float src[MAX_SAMPLE_NUM];

  for (int y=-inRadius; y<=inRadius; ++y)
  {
    for (int x=-inRadius; x<=inRadius; ++x)
    {
      float dis2 = (float)(x*x + y*y);
      if(dis2 > inRadius2)
        continue;

      const float cur_gray = out_data [(y+outRadius)*out_step +x+outRadius];
      const unsigned char cur_flag = flag_data[(y+outRadius)*flag_step+x+outRadius];
      if (cur_flag == 0)
        continue;

      const float nDirX = static_cast<float>(x);
      const float nDirY = static_cast<float>(y);
      float nOri = atan2(nDirY, nDirX);
      if (fabs(nOri - M_PI) < FLT_EPSILON)	//[-M_PI, M_PI)
      {
        nOri = static_cast<float>(-M_PI);
      }

      bool isInBound = true;
      for (int k=0; k<LIOP_NUM; k++)
      {
        const float deltaX = lsRadius * cos(nOri+k*theta);
        const float deltaY = lsRadius * sin(nOri+k*theta);

        const float sampleX = x+deltaX+outRadius;
        const float sampleY = y+deltaY+outRadius;
        float gray;

        if (!BilinearInterpolation_BorderCheck(gray, sampleX, sampleY, outPatch, flagPatch))
        {
          isInBound = false;
          break;
        }
        src[k] = gray;
      }

      if (!isInBound)
      {
        continue;
      }

      int key = 0;
      SortGray(dst, idx, src, LIOP_NUM);
      for (int k=0; k<LIOP_NUM; ++k)
      {
        key += (idx[k]+1)* m_LiopPosWeight[LIOP_NUM-k-1];
      }
      std::map<int, unsigned char>::const_iterator iter = m_LiopPatternMap.find(key);
      if (iter != m_LiopPatternMap.end())
      {
        Pixel pix;
        pix.x = x;
        pix.y = y;
        pix.f_gray = cur_gray;
        pix.weight = 1;
        pix.l_pattern = iter->second;
        pixel[pixelCount++] = pix;
      }
    }
  }

  std::sort(pixel, pixel+pixelCount, fGrayComp);	//sort by gray

  const int l_patternWidth = LIOP_NUM == 3 ? 6 : 24;
  const int dim = l_patternWidth*REGION_NUM;

  if (pixelCount >= REGION_NUM)
  {
    int curId = 0;
    int lastId = 0;
    for (int i=0; i<REGION_NUM; ++i)
    {
      const int fenceId = pixelCount*(i+1)/REGION_NUM-1;
      const float fenceGray = pixel[fenceId].f_gray;
      const int regionId = i;
      curId = lastId;

      while (true)
      {
        if (fabs(pixel[curId].f_gray-fenceGray) < FLT_EPSILON)
        {
          lastId = curId;
          break;
        }

        const int id = regionId*l_patternWidth+pixel[curId].l_pattern;
        desc[id] += pixel[curId].weight;
        ++curId;
      }

      while (true)
      {
        if(curId==pixelCount || pixel[curId].f_gray>fenceGray)
          break;

        const int id = regionId*l_patternWidth+pixel[curId].l_pattern;
        desc[id] += pixel[curId].weight;
        ++curId;
      }
    }

    ThreshNorm(desc,dim,1.f);
  }
}

void Liop_Descriptor_Extractor::extract(
  const image::Image<unsigned char> & I,
  const SIOPointFeature & feat,
  float desc[144])
{
  memset(desc, 0, sizeof(float)*144);

  //a. extract the local patch
  const int scalePatchWidth = 31;

  const int outPatchWidth = scalePatchWidth+6;
  const int outRadius = outPatchWidth/2;
  const int outRadius2 = outRadius*outRadius;
  const float scale = feat.scale();

  image::Image<float> outPatch(outPatchWidth,outPatchWidth, true, 0);
  image::Image<unsigned char> flagPatch(outPatchWidth, outPatchWidth, true, 0);

  const image::Sampler2d<image::SamplerLinear> sampler;

  // pointer alias
  unsigned char * flagPatch_data = flagPatch.data();
  float * outPatch_data = outPatch.data();

  for(int y=-outRadius; y<=outRadius; ++y)
  {
    const float ys = y*scale+feat.y();
    if(ys<0 || ys>I.Height()-1)
      continue;

    for(int x=-outRadius; x<=outRadius; ++x)
    {
      const float dis2 = (float)(x*x + y*y);
      if(dis2 > outRadius2)
        continue;

      const float xs = x*scale+feat.x();
      if(xs<0 || xs>I.Width()-1)
        continue;

      outPatch_data[(y+outRadius)*outPatchWidth+x+outRadius] = sampler(I, ys, xs);
      flagPatch_data[(y+outRadius)*outPatchWidth+x+outRadius] = 1;
    }
  }
  image::ImageGaussianFilter(image::Image<float>(outPatch), 1.2, outPatch);

  //b. creation of the LIOP ordering
  const int inRadius = scalePatchWidth/2;
  CreateLIOP_GOrder(outPatch, flagPatch, inRadius, desc);
}

template<typename T>
bool NextPermutation(std::vector<T> & p, int n)
{
  int last = n - 1;
  int i, j, k;

  i = last;
  while (i > 0 && p[i] < p[i - 1])
    i--;

  if (i == 0)
    return false;

  k = i;
  for (j = last; j >= i; j--)
    if (p[j] > p[i - 1] && p[j] < p[k])
      k = j;
  std::swap(p[k], p[i - 1]);
  for (j = last, k = i; j > k; j--, k++)
    std::swap(p[j], p[k]);

  return true;
}

void Liop_Descriptor_Extractor::GeneratePatternMap(
  std::map<int,unsigned char> & pattern_map,
  std::vector<int> & pos_weight,
  unsigned char n)
{
  pattern_map.clear();
  pos_weight.resize(n);

  //do the job
  std::vector<unsigned char> p(n);
  pos_weight[0] = 1;
  for (unsigned char i=1; i<n;i++)
  {
    pos_weight[i] = pos_weight[i-1]*10;
  }

  unsigned char count = 0;
  int key = 0;
  for (int i = 0; i < n; ++i)
  {
    p[i] = i + 1;
    key += p[i]*pos_weight[n-i-1];
  }
  pattern_map.insert(std::make_pair(key,count));
  ++count;

  while(NextPermutation(p, n))
  {
    key = 0;
    for (unsigned char i = 0; i < n; ++i)
    {
      key += p[i]*pos_weight[n-i-1];
    }
    pattern_map.insert(std::make_pair(key, count));
    ++count;
  }
}

} // namespace LIOP
} // namespace features
} // namespace openMVG
