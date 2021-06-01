//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 09:04:21 -03 2021
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
#ifndef trifocal_app_h_
#define trifocal_app_h_


struct TrifocalSampleApp {
  void ProcessCmdLine(int argc, char **argv);
  void ExtractKeypoints(); 
  void MatchKeypoints(); 
  void ComputeTracks(); 
  void Stats();
  void ExtractXYOrientation(); 
  void Display(); 
  void DisplayDesiredIds(); 
  void RobustSolve();
  void DisplayInliers();
  void DisplayInliersCamerasAndPoints(); // display inliers and and tracks
  void DisplayInliersCamerasAndPointsSIFT();
  
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
  array<Mat, io::pp::nviews> datum_; // x,y,orientation across 3 views
  array<Mat, io::pp::nviews> pxdatum_; // pixel x,y,orientation across 3 views
  // datum_[view][4 /*xy tgtx tgty*/][npts /* total number of tracks */];
  // datum_[v][1][p] = y coordinate of point p in view v
  
  // Matches
  matching::PairWiseMatches pairwise_matches_;
 
  // Tracks 
  openMVG::tracks::STLMAPTracks tracks_;
  
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers_;
};

#endif
