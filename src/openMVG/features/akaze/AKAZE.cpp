// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/akaze/AKAZE.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_filtering.hpp"
#include "openMVG/image/image_diffusion.hpp"
#include "openMVG/image/image_resampling.hpp"

#include <cmath>

namespace openMVG {
namespace features {

using namespace openMVG::image;

/// Lookup table for 2d gaussian (sigma = 2.5) where (0,0) is top left and (6,6) is bottom right
const float gauss25[7][7] = {
  {0.02546481f,  0.02350698f,  0.01849125f,  0.01239505f,  0.00708017f,  0.00344629f,  0.00142946f},
  {0.02350698f,  0.02169968f,  0.01706957f,  0.01144208f,  0.00653582f,  0.00318132f,  0.00131956f},
  {0.01849125f,  0.01706957f,  0.01342740f,  0.00900066f,  0.00514126f,  0.00250252f,  0.00103800f},
  {0.01239505f,  0.01144208f,  0.00900066f,  0.00603332f,  0.00344629f,  0.00167749f,  0.00069579f},
  {0.00708017f,  0.00653582f,  0.00514126f,  0.00344629f,  0.00196855f,  0.00095820f,  0.00039744f},
  {0.00344629f,  0.00318132f,  0.00250252f,  0.00167749f,  0.00095820f,  0.00046640f,  0.00019346f},
  {0.00142946f,  0.00131956f,  0.00103800f,  0.00069579f,  0.00039744f,  0.00019346f,  0.00008024f}
};

// Compute slice scale
static inline float Sigma( const float sigma0 , const int p , const int q , const int Q )
{
  if (p == 0 && q == 0 )
    return sigma0;
  else
    return sigma0 * powf( 2.f , p + static_cast<float>( q ) / static_cast<float>( Q ) );
}

float AKAZE::ComputeAutomaticContrastFactor( const Image<float> & src , const float percentile )
{
  const size_t nb_bin = 300;
  const int height = src.Height();
  const int width = src.Width();

  // Smooth the image
  Image<float> smoothed;
  ImageGaussianFilter( src , 1.f , smoothed , 0, 0);

  // Compute gradient
  Image<float> Lx, Ly;
  ImageScharrXDerivative( smoothed , Lx , false );
  ImageScharrYDerivative( smoothed , Ly , false );

  Image<float> & grad = smoothed; // reuse smoothed to avoid new allocation
  // grad = sqrt(Lx^2 + Ly^2)
  grad.array() = (Lx.array().square() + Ly.array().square()).sqrt();
  const float grad_max = grad.maxCoeff();

  // Compute histogram
  std::vector< size_t > histo( nb_bin, 0 );

  int nb_value = 0;

  for (int i = 1; i < height - 1; ++i )
  {
    for (int j = 1; j < width - 1; ++j )
    {
      const float val = grad( i , j );

      if (val > 0 )
      {
        int bin_id = floor( (val / grad_max ) * static_cast<float>(nb_bin) );

        // Handle overflow (need to do it in a cleaner way)
        if (bin_id == nb_bin )
          --bin_id;

        // Accumulate
        ++histo[ bin_id ];
        ++nb_value;
      }
    }
  }

  const size_t search_id = percentile * static_cast<float>(nb_value);

  size_t id_bin = 0;
  size_t acc = 0;
  while (acc < search_id && id_bin < nb_bin)
  {
    acc  += histo[ id_bin ];
    ++id_bin;
  }

  // Handle 0 bin search
  if (acc < search_id )
  {
    return 0.03f; // Only empiric value
  }
  else
  {
    return grad_max * static_cast<float>( id_bin ) / static_cast<float>( nb_bin );
  }
}

const float fderivative_factor = 1.5f;      // Factor for the multiscale derivatives

void AKAZE::ComputeAKAZESlice( const Image<float> & src , const int p , const int q , const int nbSlice ,
                        const float sigma0 , // first octave initial scale
                        const float contrast_factor ,
                        Image<float> & Li , // Diffusion image
                        Image<float> & Lx , // X derivatives
                        Image<float> & Ly , // Y derivatives
                        Image<float> & Lhess ) // Det(Hessian)
{
  const float sigma_cur = Sigma( sigma0 , p , q , nbSlice );
  const float ratio = 1 << p; //pow(2,p);
  const int sigma_scale = std::round(sigma_cur * fderivative_factor / ratio);

  Image<float> smoothed;
  if (p == 0 && q == 0 )
  {
    // Compute new image
    ImageGaussianFilter( src , sigma0 , Li, 0, 0);
  }
  else
  {
    // general case
    Image<float> in;
    if (q == 0 )  {
      ImageHalfSample( src , in );
    }
    else {
      in = src;
    }

    const float sigma_prev = ( q == 0 ) ? Sigma( sigma0 , p - 1 , nbSlice - 1 , nbSlice ) : Sigma( sigma0 , p , q - 1 , nbSlice );

    // Compute non linear timing between two consecutive slices
    const float t_prev = 0.5f * ( sigma_prev * sigma_prev );
    const float t_cur  = 0.5f * ( sigma_cur * sigma_cur );
    const float total_cycle_time = t_cur - t_prev;

    // Compute first derivatives (Scharr scale 1, non normalized) for diffusion coef
    ImageGaussianFilter( in , 1.f , smoothed, 0, 0 );

    ImageScharrXDerivative( smoothed , Lx , false );
    ImageScharrYDerivative( smoothed , Ly , false );

    // Compute diffusion coefficient
    Image<float> & diff = smoothed; // diffusivity image (reuse existing memory)
    ImagePeronaMalikG2DiffusionCoef( Lx , Ly , contrast_factor , diff );

    // Compute FED cycles
    std::vector< float > tau;
    FEDCycleTimings( total_cycle_time , 0.25f , tau );
    ImageFEDCycle( in , diff , tau );
    Li = in; // evolution image
  }

  // Compute Hessian response
  if (p == 0 && q == 0 )
  {
    smoothed = Li;
  }
  else
  {
    // Add a little smooth to image (for robustness of Scharr derivatives)
    ImageGaussianFilter( Li , 1.f , smoothed, 0, 0 );
  }

  // Compute true first derivatives
  ImageScaledScharrXDerivative( smoothed , Lx , sigma_scale );
  ImageScaledScharrYDerivative( smoothed , Ly , sigma_scale );

  // Second order spatial derivatives
  Image<float> Lxx, Lyy, Lxy;
  ImageScaledScharrXDerivative( Lx , Lxx , sigma_scale );
  ImageScaledScharrYDerivative( Lx , Lxy , sigma_scale );
  ImageScaledScharrYDerivative( Ly , Lyy , sigma_scale );

  Lx *= static_cast<float>( sigma_scale );
  Ly *= static_cast<float>( sigma_scale );

  // Compute Determinant of the Hessian
  Lhess.resize(Li.Width(), Li.Height());
  const float sigma_size_quad = Square(sigma_scale) * Square(sigma_scale);
  Lhess.array() = (Lxx.array()*Lyy.array()-Lxy.array().square())*sigma_size_quad;
}

template <typename Image>
void convert_scale(Image &src)
{
   typename Image::Tpixel min_val = src.minCoeff(), max_val = src.maxCoeff();
   src = src.array() - min_val;
   src /= max_val;
}

/// Constructor with input arguments
AKAZE::AKAZE
(
  const Image<unsigned char> & in,
  const AKAZE::Params & options
):options_(options)
{
  in_ = in.GetMat().cast<float>() / 255.f;
  options_.fDesc_factor = std::max(6.f*sqrtf(2.f), options_.fDesc_factor);
  //-- Safety check to limit the computable octave count
  const int nbOctaveMax = ceil(std::log2( std::min(in_.Width(), in_.Height())));
  options_.iNbOctave = std::min(options_.iNbOctave, nbOctaveMax);
}

/// Compute the AKAZE non linear diffusion scale space per slice
void AKAZE::Compute_AKAZEScaleSpace()
{
  float contrast_factor = ComputeAutomaticContrastFactor( in_, 0.7f );
  Image<float> input = in_;

  // Octave computation
  for (int p = 0; p < options_.iNbOctave; ++p )
  {
    contrast_factor *= (p == 0) ? 1.f : 0.75f;

    for (int q = 0; q < options_.iNbSlicePerOctave; ++q )
    {
      evolution_.emplace_back(TEvolution());
      TEvolution & evo = evolution_.back();
      // Compute Slice at (p,q) index
      ComputeAKAZESlice( input , p , q , options_.iNbSlicePerOctave , options_.fSigma0 , contrast_factor,
        evo.cur , evo.Lx , evo.Ly , evo.Lhess );

      // Prepare inputs for next slice
      input = evo.cur;

      // DEBUG octave image
#if DEBUG_OCTAVE
      std::stringstream str;
      str << "./" << "_oct_" << p << "_" << q << ".png";
      Image<float> tmp = evo.cur;
      convert_scale(tmp);
      Image< unsigned char > tmp2 ((tmp*255).cast<unsigned char>());
      WriteImage( str.str().c_str() , tmp2 );
#endif // DEBUG_OCTAVE
    }
  }
}

void detectDuplicates(
  std::vector<std::pair<AKAZEKeypoint, bool> > & previous,
  std::vector<std::pair<AKAZEKeypoint, bool> > & current)
{
  // mark duplicates - using a full search algorithm
  for (std::pair<AKAZEKeypoint, bool> & p1 : previous)
  {
    for (std::pair<AKAZEKeypoint, bool> & p2 : current)
    {
      if (p2.second == true) continue;

      // Check spatial distance
      const float dist = Square(p1.first.x - p2.first.x) + Square(p1.first.y - p2.first.y);
      if (dist <= Square(p1.first.size) && dist != 0.f)
      {
        if (p1.first.response < p2.first.response)
          p1.second = true; // mark as duplicate key point
        else
          p2.second = true; // mark as duplicate key point
        break; // no other point can be so close, so skip to the next iteration
      }
    }
  }
}

void AKAZE::Feature_Detection(std::vector<AKAZEKeypoint>& kpts) const
{
  std::vector< std::vector< std::pair<AKAZEKeypoint, bool> > > vec_kpts_perSlice(options_.iNbOctave*options_.iNbSlicePerOctave);

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int p = 0; p < options_.iNbOctave; ++p )
  {
    const float ratio = (float) (1 << p);

    for (int q = 0; q < options_.iNbSlicePerOctave; ++q )
    {
      const float sigma_cur = Sigma( options_.fSigma0 , p , q , options_.iNbSlicePerOctave );
      const Image<float> & LDetHess = evolution_[options_.iNbSlicePerOctave * p + q].Lhess;

      // Check that the point is under the image limits for the descriptor computation
      const float borderLimit =
        std::round(options_.fDesc_factor*sigma_cur*fderivative_factor/ratio)+1;

      for (int jx = borderLimit; jx < LDetHess.Height()-borderLimit; ++jx)
      for (int ix = borderLimit; ix < LDetHess.Width()-borderLimit; ++ix) {

        const float value = LDetHess(jx, ix);

        // Filter the points with the detector threshold
        if (value > options_.fThreshold &&
          value > LDetHess(jx-1, ix) &&
          value > LDetHess(jx-1, ix+1) &&
          value > LDetHess(jx-1, ix-1) &&
          value > LDetHess(jx  , ix-1) &&
          value > LDetHess(jx  , ix+1) &&
          value > LDetHess(jx+1, ix-1) &&
          value > LDetHess(jx+1, ix) &&
          value > LDetHess(jx+1, ix+1))
        {
          AKAZEKeypoint point;
          point.size = sigma_cur * fderivative_factor;
          point.octave = p;
          point.response = std::abs(value);
          point.x = ix * ratio + 0.5 * (ratio-1);
          point.y = jx * ratio + 0.5 * (ratio-1);
          point.angle = 0.0f;
          point.class_id = p * options_.iNbSlicePerOctave + q;
          vec_kpts_perSlice[options_.iNbOctave * p + q].emplace_back( point,false );
        }
      }
    }
  }

  //-- Filter duplicates
  detectDuplicates(vec_kpts_perSlice[0], vec_kpts_perSlice[0]);
  for (size_t k = 1; k < vec_kpts_perSlice.size(); ++k)
  {
    detectDuplicates(vec_kpts_perSlice[k], vec_kpts_perSlice[k]);    // detect inter scale duplicates
    detectDuplicates(vec_kpts_perSlice[k-1], vec_kpts_perSlice[k]);  // detect duplicates using previous octave
  }

  // Keep only the one marked as not duplicated
  for (size_t k = 0; k < vec_kpts_perSlice.size(); ++k)
  {
    const std::vector< std::pair<AKAZEKeypoint, bool> > & vec_kp = vec_kpts_perSlice[k];
    for (size_t i = 0; i < vec_kp.size(); ++i)
      if (!vec_kp[i].second)
        kpts.emplace_back(vec_kp[i].first);
  }
}

/// This method performs sub pixel refinement of a keypoint
bool AKAZE::Do_Subpixel_Refinement( AKAZEKeypoint & kpt, const Image<float> & Ldet) const
{
  const unsigned int ratio = (1 << kpt.octave);
  const int x = std::round(kpt.x/ratio);
  const int y = std::round(kpt.y/ratio);

  // Compute the gradient
  const float Dx = 0.5f * (Ldet(y,x+1)  - Ldet(y,x-1));
  const float Dy = 0.5f * (Ldet(y+1, x) - Ldet(y-1, x));

  // Compute the Hessian
  const float Dxx = Ldet(y, x+1) + Ldet(y, x-1) - 2.0f * Ldet(y, x);
  const float Dyy = Ldet(y+1, x) + Ldet(y-1, x) -2.0f * Ldet(y, x);
  const float Dxy = 0.25f * (Ldet(y+1, x+1) + Ldet(y-1, x-1)) - 0.25f * (Ldet(y-1, x+1) + Ldet(y+1, x-1));

  // Solve the linear system
  Eigen::Matrix<double, 2, 2> A;
  Vec2 b;
  A << Dxx, Dxy, Dxy, Dyy;
  b << -Dx, -Dy;

  const Vec2 dst = A.fullPivLu().solve(b);

  if (std::abs(dst(0)) <= 1.0 && std::abs(dst(1)) <= 1.0) {
    kpt.x += dst(0) * ratio + 0.5 * (ratio-1);
    kpt.y += dst(1) * ratio + 0.5 * (ratio-1);
    return true;
  }
  // Delete the point since its not stable
  return false;
}

/// Sub pixel refinement of the detected keypoints
void AKAZE::Do_Subpixel_Refinement(std::vector<AKAZEKeypoint>& kpts) const
{
  std::vector<AKAZEKeypoint> kpts_cpy;
  kpts_cpy.swap(kpts);
  kpts.reserve(kpts_cpy.size());

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < static_cast<int>(kpts_cpy.size()); ++i)
  {
    AKAZEKeypoint & pt = kpts_cpy[i];
    if (Do_Subpixel_Refinement(pt, this->evolution_[pt.class_id].Lhess))
    {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
      kpts.emplace_back(pt);
    }
  }
}

/// This function computes the angle from the vector given by (X Y). From 0 to 2*Pi
inline float get_angle(float x, float y)
{
  const float angle = atan2(y,x);
  // output angle between 0 and 2Pi
  return angle > 0.0f ? angle : 2.f*M_PI + angle;
}


/**
 * @brief This method computes the main orientation for a given keypoint
 * @param kpt Input keypoint
 * @note The orientation is computed using a similar approach as described in the
 * original SURF method. See Bay et al., Speeded Up Robust Features, ECCV 2006
*/
void AKAZE::Compute_Main_Orientation(
  AKAZEKeypoint& kpt,
  const Image<float> & Lx,
  const Image<float> & Ly) const
{
  int ix = 0, iy = 0, idx = 0;
  const int TABSIZE = 109;
  float resX[TABSIZE], resY[TABSIZE], Ang[TABSIZE];
  const short id[] = {6,5,4,3,2,1,0,1,2,3,4,5,6};

  // Variables for computing the dominant direction
  float sumX = 0.0f, sumY = 0.0f, max = 0.0f, ang1 = 0.0f, ang2 = 0.0f;

  // Get the information from the keypoint
  const unsigned int ratio = (1 << kpt.octave);
  const int s = std::round(kpt.size/ratio);
  const float xf = kpt.x/ratio;
  const float yf = kpt.y/ratio;

  // Calculate derivatives responses for points within radius of 6*scale
  for (int i = -6; i <= 6; ++i) {
    for (int j = -6; j <= 6; ++j) {
      if (i*i + j*j < 36) {
        iy = std::round(yf + j * s);
        ix = std::round(xf + i * s);

        const float gweight = gauss25[id[i+6]][id[j+6]];
        resX[idx] = gweight * Lx(iy, ix);
        resY[idx] = gweight * Ly(iy, ix);

        Ang[idx] = get_angle(resX[idx],resY[idx]);
        ++idx;
      }
    }
  }

  // Loop slides pi/3 window around feature point
  for (ang1 = 0.f; ang1 < 2.0f*M_PI;  ang1+=0.15f) {
    ang2 =(ang1 + M_PI / 3.0f > 2.0f * M_PI ?
      ang1 - 5.0f * M_PI / 3.0f :
      ang1 + M_PI / 3.0f);
    sumX = sumY = 0.f;

    for (int k = 0; k < idx; ++k) {
      // Get angle from the x-axis of the sample point
      const float & ang = Ang[k];

      // Determine whether the point is within the window
      if (ang1 < ang2 && ang1 < ang && ang < ang2) {
        sumX += resX[k];
        sumY += resY[k];
      }
      else if (ang2 < ang1 &&
               ((ang > 0 && ang < ang2) || (ang > ang1 && ang < 2.0f*M_PI) )) {
        sumX += resX[k];
        sumY += resY[k];
      }
    }

    // if the vector produced from this window is longer than all
    // previous vectors then this forms the new dominant direction
    if (sumX * sumX + sumY * sumY > max) {
      // store largest orientation
      max = sumX * sumX + sumY * sumY;
      kpt.angle = get_angle(sumX, sumY);
    }
  }
}

} // namespace features
} // namespace openMVG
