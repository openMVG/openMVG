// \author Ricardo FABBRI
// \date Tue Jun  1 09:04:21 -03 2021
// \author Gabriel ANDRADE \author Pierre MOULON
#ifndef trifocal_app_h_
#define trifocal_app_h_


struct TrifocalSampleApp {
 public:
   static int iteration_global_debug = 0;

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
      // TODO(gabriel): keep only 3 true tracks
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
      datum_[2].resize(4, tracks_.size()); // XXX repeat for pixdatum
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
                              cos(feature_j.orientation()), sin(feature_j.orientation());
        datum_[2].col(idx) << feature_k.x(), feature_k.y(), 
                              cos(feature_k.orientation()), sin(feature_k.orientation());
        // XXX fill up pxdatum directly WITHOUT invert_intrinsics
        for (unsigned v=0; v < 3; ++v) {
          invert_intrinsics(K_, datum_[v].col(idx).data(), datum_[v].col(idx).data()); // XXX keep datum, new std::vector: px
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

    unsigned track_id=0;
    for (const auto &track_it: tracks_) {
    //TODO: find examples of features: point in curve(3), edge(33) 
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
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawCircle(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
        svg::svgStyle().stroke("yellow", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroke("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroke("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("blue", 1));
      track_id++;
    }
    ofstream svg_file( "trifocal_track_demo.svg" );
    if (svg_file.is_open()) {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
  
  void DisplayDesiredIds() {
    //
    // Display desired ids
    //
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());

    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    unsigned track_id=0;
    for (const auto &track_it: tracks_)
    {
      bool found=false;
      for (unsigned i=0; i < n_ids; ++i)
        if (track_id == desired_ids[i])
          found = true;
          
      if (!found) {
        //cout<< found << endl;
        track_id++;
        continue;
      }
    //TODO: find examples of features: point in curve(3), edge(33) 
      auto iter = track_it.second.cbegin();
      
   uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_k = sio_regions_[2]->Features()[k];

      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroke("navy", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_k.scale(),
        svg::svgStyle().stroke("navy", 1));
      svg_stream.drawCircle(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_j.scale(),
        svg::svgStyle().stroke("navy", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroke("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroke("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("blue", 1));
      track_id++;
    }
    ofstream svg_file( "trifocal_track_desired_ids.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
//3 files trifocal_track,trifocal_inlier,track_inlier, return the correct matrices, pass to solver datum desired i,print feature sca scale
  void RobustSolve() {
    using TrifocalKernel = 
      ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                      Trifocal3PointPositionTangentialSolver>;
    Mat43 nrmdatum_; 
    Mat42 nrmdatum_;  // XXX pxdatum
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    // example: vec_inliers_ = {2, 4}  --> {33, 53} ids into orig
    array<Mat,3> Ds;
    Ds[0].resize(4,n_ids);
    Ds[1].resize(4,n_ids);
    Ds[2].resize(4,n_ids);
    //std::cerr << Ds[0].cols() << "\n";
    unsigned track_id=0;
    //going to try with calling the value of datum_[0].cols()
    for(unsigned i=0;i<datum_[0].cols();i++){
      for(unsigned j=0;j<n_ids;j++){
        if(i ==  desired_ids[j]){
          //cout << i<<"\n";
          for(unsigned k=0;k<4;k++){
            Ds[0](k, track_id) = datum_[0].col(desired_ids[j])[k];
            Ds[1](k, track_id) = datum_[1].col(desired_ids[j])[k];
            Ds[2](k, track_id) = datum_[2].col(desired_ids[j])[k];
          }
          track_id++;
        }
      }
    }
    //cout <<  Ds[0] << "\n";
    const TrifocalKernel trifocal_kernel(datum_[0], datum_[1], datum_[2], pxdatum_[0], pxdatum_[1], pxdatum_[2], K_);
    //const TrifocalKernel trifocal_kernel(Ds[0], Ds[1], Ds[2]);

    const double threshold_pix = 0.01; // 5*5 Gabriel's note : changing this for see what happens
    const unsigned max_iteration =1; // testing
    const auto model = MaxConsensus(trifocal_kernel, 
        ScorerEvaluator<TrifocalKernel>(threshold_pix), &vec_inliers_,max_iteration);
    // TODO(gabriel) recontruct from inliers and best models to show as PLY
  }

  void DisplayInliers() { //Display inliers only
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    constexpr unsigned n_inlier_pp = 3;
    unsigned desired_inliers[n_inlier_pp] = {13, 23, 63};
    unsigned track_inlier=0;
    for (const auto &track_it: tracks_)
    {
      bool inlier=false;
      for (unsigned i=0; i < n_inlier_pp; ++i)
        if (track_inlier == desired_inliers[i])
          inlier = true;
          
      if (!inlier) {
        track_inlier++;
        continue;
      }
      
      auto iter = track_it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_k = sio_regions_[2]->Features()[k];
      //cout<<"cyka"<<endl; 
      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroke("green", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
        svg::svgStyle().stroke("green", 1));
      svg_stream.drawCircle(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
        svg::svgStyle().stroke("green", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_inlier));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroke("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()),
        svg::svgStyle().stroke("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("lightblue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("lightblue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("lightblue", 1));
      track_inlier++;
    }
    ofstream svg_file( "trifocal_track_inliers.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
//display inliers and and tracks
  void DisplayInliersCamerasAndPoints() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    //constexpr unsigned n_ids_test = 5;
    //unsigned desired_ids_test[n_ids] = {13, 23, 33, 43, 53};
    // constexpr unsigned n_inlier_pp = 3;
    //unsigned desired_inliers[n_inlier_pp] = {13, 23, 43};//the only inlier that matches with desired id is 13
    vector<uint32_t>desired_inliers_vector; 
    desired_inliers_vector.resize(vec_inliers_.size()) ;
    //using this for loop for get desired_inliers_vector output
    for (unsigned j = 0; j < desired_inliers_vector.size(); j++) {
        //desired_inliers_vector.at(j) = desired_ids[vec_inliers_.at(j)];
        desired_inliers_vector.at(j) = vec_inliers_.at(j);
        //cout << desired_inliers_vector.at(j) <<" " ;
      }
    //unsigned desired_inliers[n_inlier_pp] = {desired_inliers_vector.at(13), 
    //                                         desired_inliers_vector.at(23),
    //                                         desired_inliers_vector.at(63)};//these are the selected inliers from vec_inliers_. Its easier select the result got from robustsolve() than select in robustsolve()
    unsigned track_id=0;
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
      for (unsigned i=0; i < n_ids; ++i)
        if (track_id == desired_ids[i]) { //this part is literaly overwriting the inliers
          //cout<<"blyat"<<endl;
           svg_stream.drawCircle(
              feature_i.x(), feature_i.y(), feature_i.scale(),
              svg::svgStyle().stroke("yellow", 1));
           svg_stream.drawCircle(
              feature_j.x(), feature_j.y() + images_[0].Height(), feature_k.scale(),
              svg::svgStyle().stroke("yellow", 1));
           svg_stream.drawCircle(
              feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
              svg::svgStyle().stroke("yellow", 1));
            //TODO: Tangent line segments in yellow and if inlier -> in green
            svg_stream.drawText(
              feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_id));
     
            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
              svg::svgStyle().stroke("yellow", 1)); 
            svg_stream.drawLine(
              feature_j.x(), feature_j.y() + images_[0].Height(),
              feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
              svg::svgStyle().stroke("yellow", 1));
            svg_stream.drawLine(
              feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
              feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()), //it seems that this last tangent is wrong!!
              svg::svgStyle().stroke("yellow", 1));

            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_j.x(), feature_j.y() + images_[0].Height(),
              svg::svgStyle().stroke("blue", 1));
            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_j.x(), feature_j.y() + images_[0].Height(),
              svg::svgStyle().stroke("blue", 1));
            svg_stream.drawLine(
              feature_j.x(), feature_j.y() + images_[0].Height(),
              feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
              svg::svgStyle().stroke("blue", 1));
          }
      track_id++;
    }
    
    unsigned track_inlier = 0;
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
      for (unsigned i=0; i < desired_inliers_vector.size(); ++i)
        if (track_inlier == desired_inliers_vector.at(i)) {
          //cout<<"cyka"<<endl; 
          svg_stream.drawCircle(
             feature_i.x(), feature_i.y(), feature_i.scale(),
             svg::svgStyle().stroke("green", 1));
          svg_stream.drawCircle(
            feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
            svg::svgStyle().stroke("green", 1));
          svg_stream.drawCircle(
            feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
            svg::svgStyle().stroke("green", 1));
          //TODO: Tangent line segments in yellow and if inlier -> in green
          svg_stream.drawText(
            feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_inlier));
     
          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
            svg::svgStyle().stroke("yellow", 1)); 
          svg_stream.drawLine(
            feature_j.x(), feature_j.y() + images_[0].Height(),
            feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
            svg::svgStyle().stroke("yellow", 1));
          svg_stream.drawLine(
            feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
            feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()),
            svg::svgStyle().stroke("yellow", 1));

          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_j.x(), feature_j.y() + images_[0].Height(),
            svg::svgStyle().stroke("lightblue", 1));
          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_j.x(), feature_j.y() + images_[0].Height(),
            svg::svgStyle().stroke("lightblue", 1));
          svg_stream.drawLine(
            feature_j.x(), feature_j.y() + images_[0].Height(),
            feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
            svg::svgStyle().stroke("lightblue", 1));
        }
      track_inlier++;
    }
    ofstream svg_file( "trifocal_track.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }

  void DisplayInliersCamerasAndPointsSIFT() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    unsigned track_id=0;
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    constexpr unsigned n_inlier_pp = 3;
    unsigned desired_inliers[n_inlier_pp] = {13, 23, 43};
    unsigned track_inlier=0;
    for (const auto &track_it: tracks_)
    {
      bool found=false;
      bool inlier=false;
      
      auto iter = track_it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_k = sio_regions_[2]->Features()[k];
      for (unsigned i=0; i < n_ids; ++i)
        if (track_id == desired_ids[i]) { //this part is literaly overwriting the inliers
          found = true;
      //cout<<"blyat"<<endl;//using sigma instead is a gives an error in build
      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), 2*feature_i.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_k.scale(),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawCircle(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
        svg::svgStyle().stroke("yellow", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroke("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroke("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("blue", 1));
      track_id++;
          }
      if (!found) {
        track_id++;
        continue;
      }
     track_id++;
    }
    for (const auto &track_it: tracks_)
    {
      bool inlier=false;
      
      auto iter = track_it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        k = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_k = sio_regions_[2]->Features()[k];
      for (unsigned i=0; i < n_inlier_pp; ++i)
        if (track_inlier == desired_inliers[i]){
         //cout<<"cyka"<<endl; 
         svg_stream.drawCircle(
            feature_i.x(), feature_i.y(), feature_i.scale(),
            svg::svgStyle().stroke("green", 1));
          svg_stream.drawCircle(
            feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
            svg::svgStyle().stroke("green", 1));
          svg_stream.drawCircle(
            feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(), feature_k.scale(),
            svg::svgStyle().stroke("green", 1));
          //TODO: Tangent line segments in yellow and if inlier -> in green
          svg_stream.drawText(
            feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(track_inlier));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroke("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroke("yellow", 1));
      svg_stream.drawLine(
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        feature_k.x()+ 20*sin(feature_k.orientation()), feature_k.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_k.orientation()),
        svg::svgStyle().stroke("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("lightblue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroke("lightblue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_k.x(), feature_k.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroke("lightblue", 1));
          inlier = true;
        }  
      if (!inlier) {
        track_inlier++;
        continue;
      }
     track_inlier++;
    }
    ofstream svg_file( "trifocal_track_SIFT.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
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
