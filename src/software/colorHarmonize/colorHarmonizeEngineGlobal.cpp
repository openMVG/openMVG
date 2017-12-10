// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013, 2014 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "colorHarmonizeEngineGlobal.hpp"
#include "software/SfM/SfMIOHelper.hpp"

#include "openMVG/image/image_io.hpp"
//-- Feature matches
#include <openMVG/matching/indMatch.hpp>
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/stl/stl.hpp"

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/graph/graph.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

//-- Selection Methods
#include "openMVG/color_harmonization/selection_fullFrame.hpp"
#include "openMVG/color_harmonization/selection_matchedPoints.hpp"
#include "openMVG/color_harmonization/selection_VLDSegment.hpp"

//-- Color harmonization solver
#include "openMVG/color_harmonization/global_quantile_gain_offset_alignment.hpp"

#include "openMVG/system/timer.hpp"

#include "third_party/progress/progress_display.hpp"

#include <numeric>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <functional>
#include <sstream>


namespace openMVG{

using namespace lemon;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace openMVG::lInfinity;
using namespace openMVG::sfm;

using FeatureT = features::SIOPointFeature;
using featsT = std::vector< FeatureT >;

ColorHarmonizationEngineGlobal::ColorHarmonizationEngineGlobal(
  const std::string & sSfM_Data_Filename,
  const std::string & sMatchesPath,
  const std::string & sMatchesFile,
  const std::string & sOutDirectory,
  const int selectionMethod,
  const int imgRef):
  _selectionMethod( selectionMethod ),
  _imgRef( imgRef ),
  _sMatchesFile(sMatchesFile),
  _sSfM_Data_Path(sSfM_Data_Filename),
  _sMatchesPath(sMatchesPath),
  _sOutDirectory(sOutDirectory)
{
  if ( !stlplus::folder_exists( sOutDirectory ) )
  {
    stlplus::folder_create( sOutDirectory );
  }
}

ColorHarmonizationEngineGlobal::~ColorHarmonizationEngineGlobal()
{
}

void pauseProcess()
{
  unsigned char i;
  std::cout << "\nPause : type key and press enter: ";
  std::cin >> i;
}


bool ColorHarmonizationEngineGlobal::Process()
{
  const std::string vec_selectionMethod[ 3 ] = { "fullFrame", "matchedPoints", "KVLD" };
  const std::string vec_harmonizeMethod[ 1 ] = { "quantifiedGainCompensation" };
  const int harmonizeMethod = 0;

  //-------------------
  // Load data
  //-------------------

  if ( !ReadInputData() )
    return false;
  if (_map_Matches.size() == 0 )
  {
    std::cout << std::endl << "Matches file is empty" <<std:: endl;
    return false;
  }

  //-- Remove EG with poor support:

  for (matching::PairWiseMatches::iterator iter = _map_Matches.begin();
    iter != _map_Matches.end();
    ++iter)
  {
    if (iter->second.size() < 120)
    {
      _map_Matches.erase(iter);
      iter = _map_Matches.begin();
    }
  }

  {
    graph::indexedGraph putativeGraph(getPairs(_map_Matches));

    // Save the graph before cleaning:
    graph::exportToGraphvizData(
      stlplus::create_filespec(_sOutDirectory, "input_graph_poor_supportRemoved"),
      putativeGraph);
  }

  //-------------------
  // Keep the largest CC in the image graph
  //-------------------
  if (!CleanGraph())
  {
    std::cout << std::endl << "There is no largest CC in the graph" << std::endl;
    return false;
  }

  //-------------------
  //-- Color Harmonization
  //-------------------

  //Choose image reference
  if (_imgRef == -1 )
  {
    do
    {
      std::cout << "Choose your reference image:\n";
      for (size_t i = 0; i < _vec_fileNames.size(); ++i)
      {
        std::cout << "id: " << i << "\t" << _vec_fileNames[ i ] << std::endl;
      }
    }
    while (
      !( std::cin >> _imgRef )
      || _imgRef < 0
      || _imgRef >= static_cast<int>(_vec_fileNames.size()) );
  }

  //Choose selection method
  if (_selectionMethod == -1 )
  {
    std::cout << "Choose your selection method:\n"
      << "- FullFrame: 0\n"
      << "- Matched Points: 1\n"
      << "- VLD Segment: 2\n";
    while (! ( std::cin >> _selectionMethod ) || _selectionMethod < 0 || _selectionMethod > 2)
    {
      std::cout << _selectionMethod << " is not accepted.\nPlease use a valid method number.\n";
    }
  }

  //-------------------
  // Compute remaining camera node Id
  //-------------------

  std::map<size_t, size_t> map_cameraNodeToCameraIndex; // graph node Id to 0->Ncam
  std::map<size_t, size_t> map_cameraIndexTocameraNode; // 0->Ncam correspondance to graph node Id
  std::set<size_t> set_indeximage;
  for (size_t i = 0; i < _map_Matches.size(); ++i)
  {
    matching::PairWiseMatches::const_iterator iter = _map_Matches.begin();
    std::advance(iter, i);

    const size_t I = iter->first.first;
    const size_t J = iter->first.second;
    set_indeximage.insert(I);
    set_indeximage.insert(J);
  }

  for (std::set<size_t>::const_iterator iterSet = set_indeximage.begin();
    iterSet != set_indeximage.end(); ++iterSet)
  {
    map_cameraIndexTocameraNode[std::distance(set_indeximage.begin(), iterSet)] = *iterSet;
    map_cameraNodeToCameraIndex[*iterSet] = std::distance(set_indeximage.begin(), iterSet);
  }

  std::cout << "\n Remaining cameras after CC filter : \n"
    << map_cameraIndexTocameraNode.size() << " from a total of " << _vec_fileNames.size() << std::endl;

  size_t bin      = 256;
  double minvalue = 0.0;
  double maxvalue = 255.0;

  // For each edge computes the selection masks and histograms (for the RGB channels)
  std::vector<relativeColorHistogramEdge> map_relativeHistograms[3];
  map_relativeHistograms[0].resize(_map_Matches.size());
  map_relativeHistograms[1].resize(_map_Matches.size());
  map_relativeHistograms[2].resize(_map_Matches.size());

  for (size_t i = 0; i < _map_Matches.size(); ++i)
  {
    matching::PairWiseMatches::const_iterator iter = _map_Matches.begin();
    std::advance(iter, i);

    const size_t I = iter->first.first;
    const size_t J = iter->first.second;

    const std::vector<IndMatch> & vec_matchesInd = iter->second;

    //-- Edges names:
    std::pair<std::string, std::string> p_imaNames;
    p_imaNames = make_pair( _vec_fileNames[ I ], _vec_fileNames[ J ] );
    std::cout << "Current edge : "
      << stlplus::filename_part(p_imaNames.first) << "\t"
      << stlplus::filename_part(p_imaNames.second) << std::endl;

    //-- Compute the masks from the data selection:
    Image< unsigned char > maskI ( _vec_imageSize[ I ].first, _vec_imageSize[ I ].second );
    Image< unsigned char > maskJ ( _vec_imageSize[ J ].first, _vec_imageSize[ J ].second );

    switch (_selectionMethod)
    {
      enum EHistogramSelectionMethod
      {
          eHistogramHarmonizeFullFrame     = 0,
          eHistogramHarmonizeMatchedPoints = 1,
          eHistogramHarmonizeVLDSegment    = 2,
      };

      case eHistogramHarmonizeFullFrame:
      {
        color_harmonization::commonDataByPair_FullFrame  dataSelector(
          p_imaNames.first,
          p_imaNames.second);
        dataSelector.computeMask( maskI, maskJ );
      }
      break;
      case eHistogramHarmonizeMatchedPoints:
      {
        int circleSize = 10;
        color_harmonization::commonDataByPair_MatchedPoints dataSelector(
          p_imaNames.first,
          p_imaNames.second,
          vec_matchesInd,
          _map_feats[ I ],
          _map_feats[ J ],
          circleSize);
        dataSelector.computeMask( maskI, maskJ );
      }
      break;
      case eHistogramHarmonizeVLDSegment:
      {
        color_harmonization::commonDataByPair_VLDSegment dataSelector(
          p_imaNames.first,
          p_imaNames.second,
          vec_matchesInd,
          _map_feats[ I ],
          _map_feats[ J ]);

        dataSelector.computeMask( maskI, maskJ );
      }
      break;
      default:
        std::cout << "Selection method unsupported" << std::endl;
        return false;
    }

    //-- Export the masks
    bool bExportMask = false;
    if (bExportMask)
    {
      std::string sEdge = _vec_fileNames[ I ] + "_" + _vec_fileNames[ J ];
      sEdge = stlplus::create_filespec( _sOutDirectory, sEdge );
      if ( !stlplus::folder_exists( sEdge ) )
        stlplus::folder_create( sEdge );

      std::string out_filename_I = "00_mask_I.png";
      out_filename_I = stlplus::create_filespec( sEdge, out_filename_I );

      std::string out_filename_J = "00_mask_J.png";
      out_filename_J = stlplus::create_filespec( sEdge, out_filename_J );

      WriteImage( out_filename_I.c_str(), maskI );
      WriteImage( out_filename_J.c_str(), maskJ );
    }

    //-- Compute the histograms
    Image< RGBColor > imageI, imageJ;
    ReadImage( p_imaNames.first.c_str(), &imageI );
    ReadImage( p_imaNames.second.c_str(), &imageJ );

    Histogram< double > histoI( minvalue, maxvalue, bin);
    Histogram< double > histoJ( minvalue, maxvalue, bin);

    int channelIndex = 0; // RED channel
    color_harmonization::commonDataByPair::computeHisto( histoI, maskI, channelIndex, imageI );
    color_harmonization::commonDataByPair::computeHisto( histoJ, maskJ, channelIndex, imageJ );
    relativeColorHistogramEdge & edgeR = map_relativeHistograms[channelIndex][i];
    edgeR = relativeColorHistogramEdge(map_cameraNodeToCameraIndex[I], map_cameraNodeToCameraIndex[J],
      histoI.GetHist(), histoJ.GetHist());

    histoI = histoJ = Histogram< double >( minvalue, maxvalue, bin);
    channelIndex = 1; // GREEN channel
    color_harmonization::commonDataByPair::computeHisto( histoI, maskI, channelIndex, imageI );
    color_harmonization::commonDataByPair::computeHisto( histoJ, maskJ, channelIndex, imageJ );
    relativeColorHistogramEdge & edgeG = map_relativeHistograms[channelIndex][i];
    edgeG = relativeColorHistogramEdge(map_cameraNodeToCameraIndex[I], map_cameraNodeToCameraIndex[J],
      histoI.GetHist(), histoJ.GetHist());

    histoI = histoJ = Histogram< double >( minvalue, maxvalue, bin);
    channelIndex = 2; // BLUE channel
    color_harmonization::commonDataByPair::computeHisto( histoI, maskI, channelIndex, imageI );
    color_harmonization::commonDataByPair::computeHisto( histoJ, maskJ, channelIndex, imageJ );
    relativeColorHistogramEdge & edgeB = map_relativeHistograms[channelIndex][i];
    edgeB = relativeColorHistogramEdge(map_cameraNodeToCameraIndex[I], map_cameraNodeToCameraIndex[J],
      histoI.GetHist(), histoJ.GetHist());
  }

  std::cout << "\n -- \n SOLVE for color consistency with linear programming\n --" << std::endl;
  //-- Solve for the gains and offsets:
  std::vector<size_t> vec_indexToFix;
  vec_indexToFix.push_back(map_cameraNodeToCameraIndex[_imgRef]);

  using namespace openMVG::linearProgramming;

  std::vector<double> vec_solution_r(_vec_fileNames.size() * 2 + 1);
  std::vector<double> vec_solution_g(_vec_fileNames.size() * 2 + 1);
  std::vector<double> vec_solution_b(_vec_fileNames.size() * 2 + 1);

  openMVG::system::Timer timer;

  // Red channel
  {
    OSI_CLP_SolverWrapper lpSolver(vec_solution_r.size());

    ConstraintBuilder_GainOffset cstBuilder(map_relativeHistograms[0], vec_indexToFix);
    LP_Constraints_Sparse constraint;
    cstBuilder.Build(constraint);
    lpSolver.setup(constraint);
    lpSolver.solve();
    lpSolver.getSolution(vec_solution_r);
  }
  // Green channel
  {
    OSI_CLP_SolverWrapper lpSolver(vec_solution_g.size());

    ConstraintBuilder_GainOffset cstBuilder(map_relativeHistograms[1], vec_indexToFix);
    LP_Constraints_Sparse constraint;
    cstBuilder.Build(constraint);
    lpSolver.setup(constraint);
    lpSolver.solve();
    lpSolver.getSolution(vec_solution_g);
  }
  // Blue channel
  {
    OSI_CLP_SolverWrapper lpSolver(vec_solution_b.size());

    ConstraintBuilder_GainOffset cstBuilder(map_relativeHistograms[2], vec_indexToFix);
    LP_Constraints_Sparse constraint;
    cstBuilder.Build(constraint);
    lpSolver.setup(constraint);
    lpSolver.solve();
    lpSolver.getSolution(vec_solution_b);
  }

  std::cout << std::endl
    << " ColorHarmonization solving on a graph with: " << _map_Matches.size() << " edges took (s): "
    << timer.elapsed() << std::endl
    << "LInfinity fitting error: \n"
    << "- for the red channel is: " << vec_solution_r.back() << " gray level(s)" <<std::endl
    << "- for the green channel is: " << vec_solution_g.back() << " gray level(s)" << std::endl
    << "- for the blue channel is: " << vec_solution_b.back() << " gray level(s)" << std::endl;

  std::cout << "\n\nFound solution_r:\n";
  std::copy(vec_solution_r.begin(), vec_solution_r.end(), std::ostream_iterator<double>(std::cout, " "));

  std::cout << "\n\nFound solution_g:\n";
  std::copy(vec_solution_g.begin(), vec_solution_g.end(), std::ostream_iterator<double>(std::cout, " "));

  std::cout << "\n\nFound solution_b:\n";
  std::copy(vec_solution_b.begin(), vec_solution_b.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  std::cout << "\n\nThere is :\n" << set_indeximage.size() << " images to transform." << std::endl;

  //-> convert solution to gain offset and creation of the LUT per image
  C_Progress_display my_progress_bar( set_indeximage.size() );
  for (std::set<size_t>::const_iterator iterSet = set_indeximage.begin();
    iterSet != set_indeximage.end(); ++iterSet, ++my_progress_bar)
  {
    const size_t imaNum = *iterSet;
    using Vec256 = Eigen::Matrix<double, 256, 1>;
    std::vector< Vec256 > vec_map_lut(3);

    const size_t nodeIndex = std::distance(set_indeximage.begin(), iterSet);

    const  double g_r = vec_solution_r[nodeIndex*2];
    const  double offset_r = vec_solution_r[nodeIndex*2+1];
    const  double g_g = vec_solution_g[nodeIndex*2];
    const  double offset_g = vec_solution_g[nodeIndex*2+1];
    const  double g_b = vec_solution_b[nodeIndex*2];
    const double offset_b = vec_solution_b[nodeIndex*2+1];

    for (size_t k = 0; k < 256; ++k)
    {
      vec_map_lut[0][k] = clamp( k * g_r + offset_r, 0., 255. );
      vec_map_lut[1][k] = clamp( k * g_g + offset_g, 0., 255. );
      vec_map_lut[2][k] = clamp( k * g_b + offset_b, 0., 255. );
    }

    Image< RGBColor > image_c;
    ReadImage( _vec_fileNames[ imaNum ].c_str(), &image_c );

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for (int j = 0; j < image_c.Height(); ++j)
    {
      for (int i = 0; i < image_c.Width(); ++i)
      {
        image_c(j, i)[0] = clamp(vec_map_lut[0][image_c(j, i)[0]], 0., 255.);
        image_c(j, i)[1] = clamp(vec_map_lut[1][image_c(j, i)[1]], 0., 255.);
        image_c(j, i)[2] = clamp(vec_map_lut[2][image_c(j, i)[2]], 0., 255.);
      }
    }

    const std::string out_folder = stlplus::create_filespec( _sOutDirectory,
      vec_selectionMethod[ _selectionMethod ] + "_" + vec_harmonizeMethod[ harmonizeMethod ]);
    if ( !stlplus::folder_exists( out_folder ) )
      stlplus::folder_create( out_folder );
    const std::string out_filename = stlplus::create_filespec( out_folder, stlplus::filename_part(_vec_fileNames[ imaNum ]) );

    WriteImage( out_filename.c_str(), image_c );
  }
  return true;
}

bool ColorHarmonizationEngineGlobal::ReadInputData()
{
  if ( !stlplus::is_folder( _sMatchesPath) ||
      !stlplus::is_folder( _sOutDirectory) )
  {
    std::cerr << std::endl
      << "One of the required directory is not a valid directory" << std::endl;
    return false;
  }

  if ( !stlplus::is_file( _sSfM_Data_Path ))
  {
    std::cerr << std::endl
      << "Invalid input sfm_data file: (" << _sSfM_Data_Path << ")" << std::endl;
    return false;
  }
  if (!stlplus::is_file( _sMatchesFile ))
  {
    std::cerr << std::endl
      << "Invalid match file: (" << _sMatchesFile << ")"<< std::endl;
    return false;
  }

  // a. Read input scenes views
  SfM_Data sfm_data;
  if (!Load(sfm_data, _sSfM_Data_Path, ESfM_Data(VIEWS))) {
    std::cerr << std::endl
      << "The input file \""<< _sSfM_Data_Path << "\" cannot be read" << std::endl;
    return false;
  }

  // Read images names
  for (Views::const_iterator iter = sfm_data.GetViews().begin();
    iter != sfm_data.GetViews().end(); ++iter)
  {
    const View * v = iter->second.get();
    _vec_fileNames.push_back( stlplus::create_filespec(sfm_data.s_root_path, v->s_Img_path));
    _vec_imageSize.push_back( std::make_pair( v->ui_width, v->ui_height ));
  }

  // b. Read matches

  if ( !matching::Load(_map_Matches, _sMatchesFile) )
  {
    std::cerr<< "Unable to read the geometric matrix matches" << std::endl;
    return false;
  }

  // Read features:
  for ( size_t i = 0; i < _vec_fileNames.size(); ++i )
  {
    const size_t camIndex = i;
    if ( !loadFeatsFromFile(
            stlplus::create_filespec( _sMatchesPath,
                                      stlplus::basename_part( _vec_fileNames[ camIndex ] ),
                                      ".feat" ),
            _map_feats[ camIndex ] ) )
    {
      std::cerr << "Bad reading of feature files" << std::endl;
      return false;
    }
  }

  graph::indexedGraph putativeGraph(getPairs(_map_Matches));

  // Save the graph before cleaning:
  graph::exportToGraphvizData(
      stlplus::create_filespec( _sOutDirectory, "initialGraph" ),
      putativeGraph);

  return true;
}

bool ColorHarmonizationEngineGlobal::CleanGraph()
{
  // Create a graph from pairwise correspondences:
  // - keep the largest connected component.

  graph::indexedGraph putativeGraph(getPairs(_map_Matches));

  // Save the graph before cleaning:
  graph::exportToGraphvizData(
    stlplus::create_filespec(_sOutDirectory, "initialGraph"),
    putativeGraph);

  const int connectedComponentCount = lemon::countConnectedComponents(putativeGraph.g);
  std::cout << "\n"
    << "ColorHarmonizationEngineGlobal::CleanGraph() :: => connected Component cardinal: "
    << connectedComponentCount << std::endl;

  if (connectedComponentCount > 1)  // If more than one CC, keep the largest
  {
    // Search the largest CC index
    const std::map<IndexT, std::set<lemon::ListGraph::Node>> map_subgraphs =
      openMVG::graph::exportGraphToMapSubgraphs<lemon::ListGraph, IndexT>(putativeGraph.g);
    size_t count = std::numeric_limits<size_t>::min();
    std::map<IndexT, std::set<lemon::ListGraph::Node>>::const_iterator iterLargestCC = map_subgraphs.end();
    for (std::map<IndexT, std::set<lemon::ListGraph::Node>>::const_iterator iter = map_subgraphs.begin();
        iter != map_subgraphs.end(); ++iter)
    {
      if (iter->second.size() > count)  {
        count = iter->second.size();
        iterLargestCC = iter;
      }
      std::cout << "Connected component of size : " << iter->second.size() << std::endl;
    }

    //-- Remove all nodes that are not listed in the largest CC
    for (std::map<IndexT, std::set<lemon::ListGraph::Node>>::const_iterator iter = map_subgraphs.begin();
        iter != map_subgraphs.end(); ++iter)
    {
      if (iter == iterLargestCC) // Skip this CC since it's the one we want to keep
        continue;

      const std::set<lemon::ListGraph::Node> & ccSet = iter->second;
      for (std::set<lemon::ListGraph::Node>::const_iterator iter2 = ccSet.begin();
        iter2 != ccSet.end(); ++iter2)
      {
        // Remove all outgoing edges
        for (lemon::ListGraph::OutArcIt e(putativeGraph.g, *iter2); e!=INVALID; ++e)
        {
          putativeGraph.g.erase(e);
          const IndexT Idu = (*putativeGraph.node_map_id)[putativeGraph.g.target(e)];
          const IndexT Idv = (*putativeGraph.node_map_id)[putativeGraph.g.source(e)];
          matching::PairWiseMatches::iterator iterM = _map_Matches.find(std::make_pair(Idu,Idv));
          if (iterM != _map_Matches.end())
          {
            _map_Matches.erase(iterM);
          }
          else // Try to find the opposite directed edge
          {
            iterM = _map_Matches.find(std::make_pair(Idv,Idu));
            if (iterM != _map_Matches.end())
              _map_Matches.erase(iterM);
          }
        }
      }
    }
  }

  // Save the graph after cleaning:
  graph::exportToGraphvizData(
    stlplus::create_filespec(_sOutDirectory, "cleanedGraph"),
    putativeGraph);

  std::cout << "\n"
    << "Cardinal of nodes: " << lemon::countNodes(putativeGraph.g) << "\n"
    << "Cardinal of edges: " << lemon::countEdges(putativeGraph.g) << std::endl
    << std::endl;

  return true;
}

} // namespace openMVG
