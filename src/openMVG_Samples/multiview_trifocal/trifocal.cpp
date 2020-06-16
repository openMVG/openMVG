
// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
//#include "openMVG/robust_estimation/robust_estimator_Ransac.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <memory>
#include <string>

//this is a temporary include, may be removed

#include<numeric>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::robust;
using namespace std;

struct Trifocal3PointPositionTangentialSolver {
  using trifocal_model_t = std::array<Mat34, 3>;
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };

  static void Solve(
    const Mat &bearing_0,
    const Mat &bearing_1,
    const Mat &bearing_2,
    std::vector<trifocal_model_t> *trifocal_tensor)
  {
    //...
    // TODO
    
    // io::point_tangents2params_img(p_, tgt_, tgt_ids[0], tgt_ids[1], K_, params_start_target_);
    // io::point_tangents2params(pn, tn, id_tgt0, id_tgt1, params/*[static 2*M::nparams]*/);

    // minus::solve_chicago(pn, tn, id_tgt0, id_tgt1);
  }

  static double Error(
    const trifocal_model_t &T,
    const Mat &bearing_0, // x,y,tangential, ..
    const Mat &bearing_1,
    const Mat &bearing_2)
  {
    // Return the cost related to this model and those sample data point
    // TODO(gabriel)
    
    // 1) reconstruct the 3D points and orientations
    // 2) project the 3D points and orientations on all images_
    // 3) compute error 
    
    //Gabriel: This error function is based on squared euclidian distance
   // const Vec r1 = &T - &bearing0; 
   // const Vec r2 = &T - &bearing1; 
   // const Vec r3 = &T - &bearing2;
   //const double d1 = std::inner_product(r1.begin(), r1.end(), r1.begin(), r1.end(),0);
   //const double d2 = std::inner_product(r2.begin(), r2.end(), r2.begin(), r2.end(),0);
   //const double d3 = std::inner_product(r3.begin(), r3.end(), r3.begin(), r3.end(),0);
   //return (d1+d2+d3)/6;
    return 0.0;
  }

  private:
};

static void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_coords[2], 
    double normalized_coords[2])
{
  const double *px = pix_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1]-K[1][2])/K[1][1];
  nrm[0] = (px[0] - K[0][1] - K[0][2])/K[0][0];
}

static void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_tgt_coords[3], 
    double normalized_tgt_coords[3], unsigned npts)
{
  const double *tp = pix_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*tp[1])/K[0][0];
}


template<typename SolverArg,
         typename ErrorArg,
         typename ModelArg = Trifocal3PointPositionTangentialSolver::trifocal_model_t>
class ThreeViewKernel {
 public:
   using Solver = SolverArg;
   using Model = ModelArg;
   using ErrorT = ErrorArg;

  ThreeViewKernel(const Mat &x1, const Mat &x2, const Mat &x3) : x1_(x1), x2_(x2), x3_(x3) {}

  /// The minimal number of point required for the model estimation
  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  /// The number of models that the minimal solver could return.
  enum { MAX_MODELS = Solver::MAX_MODELS };

  /// Extract required sample and fit model(s) to the sample
  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const auto
      x1 = ExtractColumns(x1_, samples),
      x2 = ExtractColumns(x2_, samples),
      x3 = ExtractColumns(x3_, samples);
    Solver::Solve(x1, x2, x3, models);
  }
  /// Return the error associated to the model and sample^nth point
  double Error(uint32_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1_.col(sample), x2_.col(sample), x3_.col(sample));
  }

  /// Number of putative point
  size_t NumSamples() const {
    return static_cast<size_t>(x1_.cols());
  }

  /// Compute a model on sampled datum
  static void Solve(const Mat &x1, const Mat &x2, const Mat &x3, std::vector<Model> *models) {
    // By offering this, Kernel types can be passed to templates.
    Solver::Solve(x1, x2, x3, models);
  }
  protected:
    const Mat &x1_, &x2_, &x3_; // corresponding point of the trifical configuration
};


struct TrifocalSampleApp {
 public:
  
   void ProcessCmdLine() {
     CmdLine cmd;
     cmd.add( make_option('a', image_filenames[0], "image_a") );
     cmd.add( make_option('b', image_filenames[1], "image_b") );
     cmd.add( make_option('c', image_filenames[2], "image_c") );
     cmd.add( make_option('K', intrinsics_filename, "K matrix") );
    
     try {
       if (argc == 1) throw std::string("Invalid command line parameter.");
       cmd.process(argc, argv);
     } catch (const std::string& s) {
       std::cerr << "Usage: " << argv[0] << '\n' << std::endl;
       std::cerr << s << std::endl;
       return EXIT_FAILURE;
     }
  }
   
  void ExtractKeypoints() {
     // Call Keypoint extractor
     using namespace openMVG::features;
     std::unique_ptr<Image_describer> image_describer;
     image_describer.reset(new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));

     if (!image_describer) {
       std::cerr << "Invalid Image_describer type" << std::endl;
       return EXIT_FAILURE;
     }
     for (const int image_idx : {0,1,2})
     {
       if (ReadImage(image_filenames_[image_idx].c_str(), &images_[image_idx]))
         image_describer->Describe(images_[image_idx], regions_per_image_[image_idx]);
     }
  }

  void MatchKeypoints() {
    //--
    // Compute corresponding points {{0,1}, {1,2}}
    //--
    //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
    //std::unique_ptr<Matcher> collectionMatcher(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    //collectionMatcher->Match(regions_provider, {{0,1}, {1,2}}, pairwise_matches_, &progress);
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      * regions_per_image_.at(0).get(),
      * regions_per_image_.at(1).get(),
      pairwise_matches_[{0,1}]);
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      * regions_per_image_.at(1).get(),
      * regions_per_image_.at(2).get(),
      pairwise_matches_[{1,2}]);
  }

  void ComputeTracks() {
      openMVG::tracks::TracksBuilder track_builder;
      track_builder.Build(pairwise_matches_);
      track_builder.Filter(3);
      track_builder.ExportToSTL(tracks_);
  }

  void Stats() {
      // Display some statistics
      std::cout
        <<  regions_per_image_.at(0)->RegionCount() << " #Features on image A" << std::endl
        <<  regions_per_image_.at(1)->RegionCount() << " #Features on image B" << std::endl
        <<  regions_per_image_.at(2)->RegionCount() << " #Features on image C" << std::endl
        << pairwise_matches_.at({0,1}).size() << " #matches with Distance Ratio filter" << std::endl
        << pairwise_matches_.at({1,2}).size() << " #matches with Distance Ratio filter" << std::endl
        << tracks_.size() << " #tracks" << std::endl;
  }

  void ExtractXYOrientation() {
      sio_regions_ = std::array<const SIFT_Regions*, 3> ({
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(0).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(1).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(2).get())
      });
      //
      // Build datum (corresponding {x,y,orientation})
      //
      std::array<Mat, 3> datum;              
      datum[0].resize(4, tracks_.size());
      datum[1].resize(4, tracks_.size());
      datum[2].resize(4, tracks_.size());
      int idx = 0;
      for (const auto &track_it: tracks_)
      {
        auto iter = track_it.second.cbegin();
        const uint32_t
          i = iter->second,
          j = (++iter)->second,
          k = (++iter)->second;
        //
        const auto feature_i = sio_regions[0]->Features()[i];
        const auto feature_j = sio_regions[1]->Features()[j];
        const auto feature_k = sio_regions[2]->Features()[k];
        datum[0].col(idx) << feature_i.x(), feature_i.y(), std::cos(feature_i.orientation()), std::sin(feature_i.orientation());
        // datum[0].col(idx) << feature_i.x(), feature_i.y(), feature_i.orientation();
        datum[1].col(idx) << feature_j.x(), feature_j.y(), std::cos(feature_i.orientation()), std::sin(feature_i.orientation());
    feature_j.orientation();
        datum[2].col(idx) << feature_k.x(), feature_k.y(), std::cos(feature_i.orientation()), std::sin(feature_i.orientation());
        
        //Gabriel:Calling both K invertions:
        //
        std::cout << invert_intrinsics(K, datum[0].col(idx).data(), datum[0].col(idx).data(), tracks_.size());
        // std::cout << invert_intrinsics_tgt((double *)[3][3](K.data()), datum[0].col(idx).data()+2, datum[0].col(idx).data()+2, tracks_.size());
        ++idx;
      }
      std::cout <<  datum[0] << "blyat" << std::endl;
  }

  void Display() {
    //
    // Display demo
    //
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());

    for (const auto &track_it: tracks_)
    {
      auto iter = track_it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions[0]->Features()[i];
      const auto feature_j = sio_regions[1]->Features()[j];
      const auto feature_k = sio_regions[2]->Features()[k];

      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_k.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawCircle(
        feature_k.x(), feature_j.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_j.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("blue", 1));
    }
    std::ofstream svg_file( "trifocal_track.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }

  void RobustSolve() {
    using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, Trifocal3PointPositionTangentialSolver>;
    const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

    const double threshold_pix = 2.0;
    std::vector<uint32_t> vec_inliers;
    const auto model = MaxConsensus(trifocal_kernel, ScorerEvaluator<TrifocalKernel>(threshold_pix), &vec_inliers);
  }

  void DisplayInliersCamerasAndPoints() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

  }

  // ---------------------------------------------------------------------------
  // Data
  //
  //
  // 3x3 intrinsic matrix for this default test case
  // This representation is specific for fast non-homog action
  // Just eliminate last row 
  //
  // This matrix is calib.intrinsic for the synthcurves spherical dataset
  double K_[2][3] = {  // some default value for testing
    {2584.9325098195013197, 0, 249.77137587221417903},
    {0, 2584.7918606057692159, 278.31267937919352562}
   //  0 0 1 
  };
  
  // The three images_ used to compute the trifocal configuration
  std::array<std::string, 3> image_filenames_;
  std::array<std::string> intrinsics_filename_;
  std::array<Image<unsigned char>, 3> images_;
  
  // Features
  std::map<IndexT, std::unique_ptr<features::Regions>> regions_per_image_;
  std::array<const SIFT_Regions*, 3> sio_regions_; // a cast on regions_per_image_
  
  // Matches
  matching::PairWiseMatches pairwise_matches_;
 
  // Tracks 
  openMVG::tracks::STLMAPTracks tracks_;
}

int main(int argc, char **argv) {
  TrifocalSampleApp T;
  
  T.ProcessCmdLine();
//  T.ExtractKeypoints();
//  T.MatchKeypoints();
//  T.Stats();
//  T.ExtractXYOrientation();
//  T.Display();
//  T.RobustSolve();
//  T.DisplayInliersCamerasAndPoints();

  return EXIT_SUCCESS;
}
