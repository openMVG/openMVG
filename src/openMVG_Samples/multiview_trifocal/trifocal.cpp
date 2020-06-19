
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

//these are temporary includes, may be removed

#include<numeric>
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

// Mat is Eigen::MatrixXd - matrix of doubles with dynamic size
// Vec3 is Eigen::Vector3d - Matrix< double, 3, 1 >

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::robust;
using namespace std;
using SIFT_Regions = openMVG::features::SIFT_Regions;

struct Trifocal3PointPositionTangentialSolver {
  using trifocal_model_t = array<Mat34, 3>;
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };

  // datum_i[4 /*xy tgtx tgty*/][pp:npoints /* 3 for Chicago */]
  static void Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      vector<trifocal_model_t> *trifocal_tensor) {

    double p[pp::nviews][pp::npoints][io::ncoords2d];
    double tgt[pp::nviews][pp::npoints][io::ncoords2d]; 
    
    // pack into solver's efficient representation
    for (unsigned ip=0; ip < pp:npoints; ++ip) {
      p[0][ip][0] = datum_0(0,ip);
      p[0][ip][1] = datum_0(1,ip);
      tgt[0][ip][0] = datum_0(2,ip); 
      tgt[0][ip][1] = datum_0(3,ip); 
      
      p[1][ip][0] = datum_1(0,ip);
      p[1][ip][1] = datum_1(1,ip);
      tgt[1][ip][0] = datum_1(2,ip); 
      tgt[1][ip][1] = datum_1(3,ip); 
      
      p[2][ip][0] = datum_2(0,ip);
      p[2][ip][1] = datum_2(1,ip);
      tgt[2][ip][0] = datum_2(2,ip); 
      tgt[2][ip][1] = datum_2(3,ip); 
    }
    
    unsigned nsols_final = 0;
    unsigned id_sols[M::nsols];
    double  solutions_cams[M::nsols][pp::nviews-1][4][3];  // first camera is always [I | 0]
    
    minus<chicago>::solve(p, tgt, cameras, id_sols, &nsols_final);

    // xxx
    for (unsigned s=0; s < nsols_final; ++s)
      trifocal_tensor->push_back(cameras[id_sols[s]]);
  }
  
  // Gabriel's comment: If bearing is the bearing vector of the camera, Vec3 should be used instead of Mat32 or use &bearing.data()[0] 
  static double Error(
    const trifocal_model_t &T,
    const Vec3 &bearing_0, // x,y,tangential, ..
    const Vec3 &bearing_1,
    const Vec3 &bearing_2) {
    // Return the cost related to this model and those sample data point
    // TODO(gabriel)
    
    // 1) reconstruct the 3D points and orientations
    // 2) project the 3D points and orientations on all images_
    // 3) compute error 
    
    // Gabriel's note: I'm using triangulation.hpp file, for best usage I'll create an Vec3 array
    const array<Vec3, 3>triangulated;
    openMVG::TriangulateDLT( T[0], bearing_0, T[1], bearing_1, triangulated[0] );
    openMVG::TriangulateDLT( T[1], bearing_1, T[2], bearing_2, triangulated[1] );
    openMVG::TriangulateDLT( T[2], bearing_2, T[0], bearing_0, triangulated[2] );
    
    // Computing the projection of triangulated points using projection.hpp
    const array<Vec2, 3>projected;
    projected[0] << openMVG::Project( T[0], triangulated[0] );
    projected[1] << openMVG::Project( T[1], triangulated[1] );
    projected[2] << openMVG::Project( T[2], triangulated[2] );
    
    //computing the projection of inputs
    const array<Vec2, 3> real_projection;
    real_projection[0] << bearing_0[0]/bearing_0[2],
                          bearing_0[1]/bearing_0[2];
    real_projection[1] << bearing_1[0]/bearing_1[2], 
                          bearing_1[1]/bearing_1[2];
    real_projection[2] << bearing_2[0]/bearing_2[2], 
                          bearing_2[1]/bearing_2[2];

    //computing the euclidian distance and error
    const array<double, 3> distance;
    for ( unsigned i =0; i<3 ; i++){
       distance[i] = inner_product( (projected[i]-real_projection[i]), (projected[i]-real_projection[i])+3, (projected[i]-real_projection[i]), 0.0 );
    }
    //distance[0] = inner_product( (projected[0]-real_projection[0])[0], (projected[0]-real_projection[0])[2], (projected[0]-real_projection[0])[0], 0.0 );
    //distance[1] = inner_product( (projected[1]-real_projection[1])[0], (projected[1]-real_projection[1])[2], (projected[1]-real_projection[1])[0], 0.0 );
    //distance[2] = inner_product( (projected[2]-real_projection[2])[0], (projected[2]-real_projection[2])[2], (projected[2]-real_projection[2])[0], 0.0 );
    double sum = 0.0;
    for ( unsigned i = 0; i < 3; i++){
          sum += distance[i];
    }

    //Gabriel's comment: Using square Euclidian distance
    return sum/6;
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
  nrm[1] = (px[1] - K[1][2]) /K[1][1];
  nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
}

static void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_tgt_coords[2], 
    double normalized_tgt_coords[2])
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
  void Fit(const vector<uint32_t> &samples, vector<Model> *models) const {
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

  /// Compute a model on sampled datum_
  static void Solve(const Mat &x1, const Mat &x2, const Mat &x3, vector<Model> *models) {
    // By offering this, Kernel types can be passed to templates.
    Solver::Solve(x1, x2, x3, models);
  }
  protected:
    const Mat &x1_, &x2_, &x3_; // corresponding point of the trifical configuration
    // x_i[4 /*xy tgtx tgty*/][npts /* total number of tracks */]
};


struct TrifocalSampleApp {
 public:
  
   void ProcessCmdLine(int argc, char **argv) {
     CmdLine cmd;
     cmd.add( make_option('a', image_filenames_[0], "image_a") );
     cmd.add( make_option('b', image_filenames_[1], "image_b") );
     cmd.add( make_option('c', image_filenames_[2], "image_c") );
     cmd.add( make_option('K', intrinsics_filename_, "K matrix") );
    
     try {
       if (argc == 1) throw string("Invalid command line parameter.");
       cmd.process(argc, argv);
     } catch (const string& s) {
       cerr << "Usage: " << argv[0] << '\n' << endl;
       cerr << s << endl;
       exit(EXIT_FAILURE);
     }
  }
   
  void ExtractKeypoints() {
     // Call Keypoint extractor
     using namespace openMVG::features;
     unique_ptr<Image_describer> image_describer;
     image_describer.reset(new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));

     if (!image_describer) {
       cerr << "Invalid Image_describer type" << endl;
       exit(EXIT_FAILURE);
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
    //unique_ptr<Matcher> collectionMatcher(new Cascade_Hashing_Matcher_Regions(fDistRatio));
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
      cout
        <<  regions_per_image_.at(0)->RegionCount() << " #Features on image A" << endl
        <<  regions_per_image_.at(1)->RegionCount() << " #Features on image B" << endl
        <<  regions_per_image_.at(2)->RegionCount() << " #Features on image C" << endl
        << pairwise_matches_.at({0,1}).size() << " #matches with Distance Ratio filter" << endl
        << pairwise_matches_.at({1,2}).size() << " #matches with Distance Ratio filter" << endl
        << tracks_.size() << " #tracks" << endl;
  }

  void ExtractXYOrientation() {
      sio_regions_ = array<const SIFT_Regions*, 3> ({
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(0).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(1).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(2).get())
      });
      //
      // Build datum_ (corresponding {x,y,orientation})
      //
      datum_[0].resize(4, tracks_.size());
      datum_[1].resize(4, tracks_.size());
      datum_[2].resize(4, tracks_.size());
      int idx = 0;
      for (const auto &track_it: tracks_) {
        auto iter = track_it.second.cbegin();
        const uint32_t
          i = iter->second,
          j = (++iter)->second,
          k = (++iter)->second;
        //
        const auto feature_i = sio_regions_[0]->Features()[i];
        const auto feature_j = sio_regions_[1]->Features()[j];
        const auto feature_k = sio_regions_[2]->Features()[k];
        datum_[0].col(idx) << feature_i.x(), feature_i.y(), 
          cos(feature_i.orientation()), sin(feature_i.orientation());
        // datum_[0].col(idx) << feature_i.x(), feature_i.y(), feature_i.orientation();
        datum_[1].col(idx) << feature_j.x(), feature_j.y(), 
          cos(feature_i.orientation()), sin(feature_i.orientation());
        datum_[2].col(idx) << feature_k.x(), feature_k.y(), 
          cos(feature_i.orientation()), sin(feature_i.orientation());
        for (unsigned v=0; v < 3; ++v) {
          invert_intrinsics(K_, datum_[v].col(idx).data(), datum_[v].col(idx).data());
          invert_intrinsics_tgt(K_, datum_[v].col(idx).data()+2, datum_[v].col(idx).data()+2);
        }
        ++idx;
      }
  }

  void Display() {
    //
    // Display demo
    //
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());

    for (const auto &track_it: tracks_)
    {
      auto iter = track_it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_k = sio_regions_[2]->Features()[k];

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
    ofstream svg_file( "trifocal_track.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }

  void RobustSolve() {
    using TrifocalKernel = 
      ThreeViewKernel<Trifocal3PointPositionTangentialSolver, Trifocal3PointPositionTangentialSolver>;
    const TrifocalKernel trifocal_kernel(datum_[0], datum_[1], datum_[2]);

    const double threshold_pix = 2.0;
    const auto model = MaxConsensus(trifocal_kernel, ScorerEvaluator<TrifocalKernel>(threshold_pix), &vec_inliers_);
  }

  void DisplayInliersCamerasAndPoints() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

  }

  
  // ---------------------------------------------------------------------------
  // Data
  // ---------------------------------------------------------------------------
 

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
  
  // The three images used to compute the trifocal configuration
  array<string, 3> image_filenames_;
  string intrinsics_filename_;
  array<Image<unsigned char>, 3> images_;
  
  // Features
  map<IndexT, unique_ptr<features::Regions>> regions_per_image_;
  array<const SIFT_Regions*, 3> sio_regions_; // a cast on regions_per_image_
  array<Mat, pp:nviews> datum_; // x,y,orientation across 3 views
  // datum_[view][4 /*xy tgtx tgty*/][npts /* total number of tracks */];
  // datum_[v][1][p] = y coordinate of point p in view v
  
  // Matches
  matching::PairWiseMatches pairwise_matches_;
 
  // Tracks 
  openMVG::tracks::STLMAPTracks tracks_;
  
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers_;
};

int main(int argc, char **argv) {
  TrifocalSampleApp T;
  
  T.ProcessCmdLine(argc, argv);
  T.ExtractKeypoints();
  T.MatchKeypoints();
  T.ComputeTracks();
  T.Stats();
  T.ExtractXYOrientation();
  T.Display();
  T.RobustSolve();
  T.DisplayInliersCamerasAndPoints();

  return EXIT_SUCCESS;
}
