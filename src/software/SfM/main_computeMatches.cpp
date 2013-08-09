
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "patented/sift/SIFT.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include "software/SfM/pairwiseAdjacencyDisplay.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace std;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sImaExtension;
  std::string sOutDir = "";
  float fDistRatio = .6f;
  bool bOctMinus1 = false;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('e', sImaExtension, "ext") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', fDistRatio, "distratio") );
  cmd.add( make_option('s', bOctMinus1, "octminus1") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imadir path] "
      << "[-e|--ext extension '*.jpg' or '*.png'] "
      << "[-o|--outdir path] "
      << "[-r|--distratio 0.6] "
      << "[-s|--octminus1 0 or 1] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << "main_computeMatches" << std::endl
            << "--imadir " << sImaDirectory << std::endl
            << "--ext " << sImaExtension << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--distratio " << fDistRatio << std::endl
            << "--octminus1 " << bOctMinus1 << std::endl;

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // -----------------------------
  // a. List images
  // b. Compute features and descriptor
  // c. Compute putatives descriptor matches
  // d. Geometric filtering of putatives matches
  // e. Export some statistics
  // -----------------------------

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  //---------------------------------------
  // a. List images
  //---------------------------------------

  std::vector<std::string> vec_fileNames;
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;

  if (!stlplus::folder_exists(sImaDirectory)) {
    std::cerr << "It is an invalid input image directory" << std::endl;
    return EXIT_FAILURE;
  } else  {

    //---------------------------------------
    // Look for images in the given directory
    //---------------------------------------
    vec_fileNames = stlplus::folder_wildcard(sImaDirectory,
      sImaExtension, false, true);

    std::sort(vec_fileNames.begin(), vec_fileNames.end());

    std::ofstream file(stlplus::create_filespec(sOutDir, "/lists", ".txt").c_str());
    std::copy(vec_fileNames.begin(), vec_fileNames.end(),
              std::ostream_iterator<std::string>(file, "\n"));
    file.close();

    for (size_t i=0; i < vec_fileNames.size(); ++i)  {
      vec_fileNames[i] = stlplus::create_filespec(
        stlplus::folder_append_separator(sImaDirectory), vec_fileNames[i]);
    }
    // DEBUG INFO
    std::cout << std::endl << "IMAGE(S) :" << std::endl;
    copy(vec_fileNames.begin(), vec_fileNames.end(), ostream_iterator<string>(cout, "\n"));

    if (vec_fileNames.empty())
    {
      std::cout << "\n No images in the provided directory.";
      return EXIT_FAILURE;
    }
  }


  //---------------------------------------
  // b. Compute features and descriptor
  //    - extract sift features and descriptor
  //    - if keypoints already computed, re-load them
  //    - else save features and descriptors and disk
  //---------------------------------------

  typedef Descriptor<float, 128> descT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> featsT;
  typedef vector<descT > descsT;
  typedef KeypointSet<featsT, descsT > KeypointSetT;

  vec_imagesSize.resize(vec_fileNames.size());
  {
    // extract SIFT features
    cout << endl << endl << "EXTRACT SIFT FEATURES" << endl;
    C_Progress_display my_progress_bar( vec_fileNames.size() );
    for(size_t i=0; i < vec_fileNames.size(); ++i)  {
      KeypointSetT kpSet;
      Image<RGBColor> imageRGB;
      Image<unsigned char> imageGray;

      std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "desc");

      //Test if descriptor and feature was already computed
      if (stlplus::file_exists(sFeat) && stlplus::file_exists(sDesc)) {

        ReadImage(vec_fileNames[i].c_str(), &imageRGB);
        vec_imagesSize[i] = make_pair(imageRGB.Width(), imageRGB.Height());
      }
      else  { //Not already computed, so compute and save

        ReadImage(vec_fileNames[i].c_str(), &imageRGB);
        Rgb2Gray(imageRGB, &imageGray);
        // Compute features and descriptors and export them to file
        SIFTDetector(imageGray,  kpSet.features(), kpSet.descriptors(), bOctMinus1);
        kpSet.saveToBinFile(sFeat, sDesc);
        vec_imagesSize[i] = make_pair(imageRGB.Width(), imageRGB.Height());
      }
      ++my_progress_bar;
    }
  }

  //-- Load descriptor in memory
  std::map<size_t, KeypointSetT > map_featAndDesc; //empty descriptor
  std::map<size_t, descT::bin_type * > map_Desc; // descriptor as contiguous memory
  {
     for (size_t j = 0; j < vec_fileNames.size(); ++j)  {
        // Load descriptor of Jnth image
        std::string sFeatJ = stlplus::create_filespec(sOutDir,
          stlplus::basename_part(vec_fileNames[j]), "feat");

        std::string sDescJ = stlplus::create_filespec(sOutDir,
          stlplus::basename_part(vec_fileNames[j]), "desc");
        loadFeatsFromFile(sFeatJ, map_featAndDesc[j].features());

        KeypointSetT kpSetI;
        loadDescsFromBinFile(sDescJ, kpSetI.descriptors());
        map_Desc[j] = new descT::bin_type[kpSetI.descriptors().size() * descT::static_size];
        for(size_t k=0; k < kpSetI.descriptors().size(); ++k)
          memcpy(
                 &map_Desc[j][k*descT::static_size],
                 kpSetI.descriptors()[k].getData(),
                 descT::static_size*sizeof(descT::bin_type));
     }
  }

  //---------------------------------------
  // c. Compute putatives descriptor matches
  //    - L2 descriptor matching
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------
  typedef std::map< std::pair<size_t, size_t>, std::vector<IndMatch> > IndexedMatchPerPair;

  IndexedMatchPerPair map_PutativesMatches;
  if (stlplus::file_exists(sOutDir + "/matches.putative.txt"))
  {
    PairedIndMatchImport(sOutDir + "/matches.putative.txt",map_PutativesMatches);
    cout << endl << endl << "PUTATIVE MATCHES -- PREVIOUS RESULTS LOADED" << endl;
  }
  else
  {
    std::cout << std::endl << std::endl << "PUTATIVE MATCHES" << std::endl;
#ifdef USE_OPENMP
    std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
    C_Progress_display my_progress_bar2( vec_fileNames.size()*(vec_fileNames.size()-1) / 2.0 );

    for (size_t i = 0; i < vec_fileNames.size(); ++i)
    {
      // Define the matcher and the used metric (Squared L2)
      // ANN matcher could be defined as follow:
      typedef flann::L2<descT::bin_type> MetricT;
      typedef ArrayMatcher_Kdtree_Flann<descT::bin_type, MetricT> MatcherT;
      // Brute force matcher is defined as following:
      //typedef L2_Vectorized<descT::bin_type> MetricT;
      //typedef ArrayMatcherBruteForce<descT::bin_type, MetricT> MatcherT;

      // Load descriptor of Inth image
      const KeypointSetT & kpSetI = map_featAndDesc[i];
      const descT::bin_type * tab0 = map_Desc[i];

      MatcherT matcher10;
      ( matcher10.Build(tab0, kpSetI.features().size(), descT::static_size) );

#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
#endif
      for (int j = i+1; j < (int)vec_fileNames.size(); ++j)
      {
        // Load descriptor of Jnth image
        const KeypointSetT & kpSetJ = map_featAndDesc[j];
        const descT::bin_type * tab1 = map_Desc[j];

        const size_t NNN__ = 2;
        std::vector<int> vec_nIndice10;
        std::vector<MetricT::ResultType> vec_fDistance10;

        //Find left->right
        matcher10.SearchNeighbours(tab1, kpSetJ.features().size(), &vec_nIndice10, &vec_fDistance10, NNN__);

        std::vector<IndMatch> vec_FilteredMatches;
        std::vector<int> vec_NNRatioIndexes;
        NNdistanceRatio( vec_fDistance10.begin(), // distance start
          vec_fDistance10.end(),  // distance end
          NNN__, // Number of neighbor in iterator sequence (minimum required 2)
          vec_NNRatioIndexes, // output (index that respect Lowe Ratio)
          Square(fDistRatio)); // squared dist ratio due to usage of a squared metric

        for (size_t k=0; k < vec_NNRatioIndexes.size()-1&& vec_NNRatioIndexes.size()>0; ++k)
        {
          vec_FilteredMatches.push_back(
            IndMatch(vec_nIndice10[vec_NNRatioIndexes[k]*NNN__],
                     vec_NNRatioIndexes[k]) );
        }

        // Remove duplicates
        IndMatch::getDeduplicated(vec_FilteredMatches);

        // Remove matches that have the same X,Y coordinates
        {
          IndMatchDecorator<float> matchDeduplicator(
            vec_FilteredMatches, kpSetI.features(), kpSetJ.features());
          matchDeduplicator.getDeduplicated(vec_FilteredMatches);

#ifdef USE_OPENMP
  #pragma omp critical
#endif
          {
            map_PutativesMatches.insert( make_pair( make_pair(i,j), vec_FilteredMatches ));
          }

        }
        ++my_progress_bar2;
      }
    }
    //---------------------------------------
    //-- Export putative matches
    //---------------------------------------
    std::ofstream file (std::string(sOutDir + "/matches.putative.txt").c_str());
    if (file.is_open())
      PairedIndMatchToStream(map_PutativesMatches, file);
    file.close();
  }

  //-- Free descriptor memory:
  for (std::map<size_t, descT::bin_type * >::const_iterator itDesc = map_Desc.begin();
    itDesc != map_Desc.end(); ++itDesc)
  {
    delete [] itDesc->second;
  }
  map_Desc.clear();

  //---------------------------------------
  // d. Geometric filtering of putatives matches
  //    - AContrario Estimation of the Fundamental matrix
  //    - Use a upper bound for the plausible F matrix
  //      acontrario estimated threshold
  //---------------------------------------


  IndexedMatchPerPair map_GeometricMatches_F;
  {
    cout << endl << endl << " - GEOMETRIC FILTERING - " << endl;
    C_Progress_display my_progress_bar3( map_PutativesMatches.size() );

#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < (int)map_PutativesMatches.size(); ++i)
    {
      IndexedMatchPerPair::const_iterator iter = map_PutativesMatches.begin();
      advance(iter,i);

      size_t iIndex = iter->first.first;
      size_t jIndex = iter->first.second;
      const vector<IndMatch> & vec_FilteredMatches = iter->second;

      // Load features of Inth image
      const std::vector<FeatureT> & kpSetI = map_featAndDesc[iIndex].features();
      // Load features of Jnth image
      const std::vector<FeatureT> & kpSetJ = map_featAndDesc[jIndex].features();

      //-- Copy point to array in order to estimate fundamental matrix :
      const size_t n = vec_FilteredMatches.size();
      Mat xA(2,n), xB(2,n);

      for (size_t i=0; i < vec_FilteredMatches.size(); ++i)  {
        const FeatureT & imaA = kpSetI[vec_FilteredMatches[i]._i];
        const FeatureT & imaB = kpSetJ[vec_FilteredMatches[i]._j];
        xA.col(i) = imaA.coords().cast<double>();
        xB.col(i) = imaB.coords().cast<double>();
      }

      //-- Fundamental matrix robust estimation
      {
        std::vector<size_t> vec_inliers;
        // Define the AContrario adapted Fundamental matrix solver
        typedef ACKernelAdaptor<
          openMVG::fundamental::kernel::SevenPointSolver,
          openMVG::fundamental::kernel::SimpleError,
          UnnormalizerT,
          Mat3>
          KernelType;

        KernelType kernel(xA, vec_imagesSize[iIndex].first, vec_imagesSize[iIndex].second,
                          xB, vec_imagesSize[jIndex].first, vec_imagesSize[jIndex].second, true);

        // Robustly estimate the Fundamental matrix with A Contrario ransac
        Mat3 F;
        double upper_bound_precision = 4.0; // upper_bound of 4 pixels
        std::pair<double,double> ACRansacOut =
          ACRANSAC(kernel, vec_inliers, 4096, &F, upper_bound_precision);
        const double & threshold = ACRansacOut.first;
        const double & NFA = ACRansacOut.second;

        if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)  {
          std::vector<IndMatch> vec_matchesF;
          vec_matchesF.reserve(vec_inliers.size());
          for (size_t i=0; i < vec_inliers.size(); ++i)  {
            vec_matchesF.push_back( vec_FilteredMatches[vec_inliers[i]] );
          }
#ifdef USE_OPENMP
  #pragma omp critical
#endif
          {
            map_GeometricMatches_F[make_pair(iIndex,jIndex)] = vec_matchesF;
          }
        }
      }
      ++my_progress_bar3;
    }
    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    std::ofstream file (string(sOutDir + "/matches.f.txt").c_str());
    if (file.is_open())
      PairedIndMatchToStream(map_GeometricMatches_F, file);
    file.close();

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's Epipolar matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
      map_GeometricMatches_F,
      stlplus::create_filespec(sOutDir, "EpipolarAdjacencyMatrix", "svg"));
  }
  return EXIT_SUCCESS;
}


