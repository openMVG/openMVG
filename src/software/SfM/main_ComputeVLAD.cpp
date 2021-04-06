// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/graph.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching_image_collection/Retrieval_Helpers.hpp"
#include "openMVG/matching_image_collection/Vlad.hpp"
#include "openMVG/sfm/pipelines/sfm_preemptive_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::retrieval;

using namespace std;
using namespace svg;

// TODO: these two output function could be factored and put somewhere else
// Display the image retrieval matrix
template <typename Order>
void saveRetrievalMatrix(
    const string &filename, const SfM_Data &sfm_data,
    const IndexedPairwiseSimilarity<Order> &result_ordered_by_similarity) {
  const size_t num_neighbors =
      result_ordered_by_similarity.begin()->second.size();

  const int size = 1024;
  svgDrawer svg_stream(size * num_neighbors + 1,
                       size * result_ordered_by_similarity.size());
  int x_offset = 0;
  int y_offset = 0;

  for (const auto result_it : result_ordered_by_similarity) {
    const auto ref_view_id = result_it.first;
    const auto ref_view = sfm_data.GetViews().at(ref_view_id);
    const std::string ref_view_filename =
        stlplus::create_filespec(sfm_data.s_root_path, ref_view->s_Img_path);

    svg_stream.drawImage(ref_view_filename, size, size, x_offset * size,
                         y_offset * size);

    const auto &retrieval_list = result_it.second;
    for (const auto retrieval_list_it : retrieval_list) {
      ++x_offset;
      const auto found_view_id = retrieval_list_it.second;

      const auto found_view = sfm_data.GetViews().at(found_view_id);
      const std::string found_view_filename = stlplus::create_filespec(
          sfm_data.s_root_path, found_view->s_Img_path);

      svg_stream.drawImage(found_view_filename, size, size, x_offset * size,
                           y_offset * size);
    }
    x_offset = 0;
    ++y_offset;
  }

  ofstream svg_file(filename.c_str());
  svg_file << svg_stream.closeSvgFile().str();
  svg_file.close();
}

void saveAdjacencyMatrixViewGraph(const Pair_Set &resulting_pairs,
                                  const SfM_Data &sfm_data,
                                  const std::string &base_folder) {
  PairWiseMatches vlad_pairwise_matches;
  for (const auto pair_it : resulting_pairs) {
    vlad_pairwise_matches[pair_it] = {
        0};  // Set at last a match to be displayed
  }
  PairWiseMatchingToAdjacencyMatrixSVG(
      sfm_data.GetViews().size(), vlad_pairwise_matches,
      stlplus::create_filespec(base_folder, "PutativeAdjacencyMatrixRetrieval", "svg"));

  //-- export view pair graph once putative graph matches have been computed
  {
    std::set<IndexT> set_viewIds;
    std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
                   std::inserter(set_viewIds, set_viewIds.begin()),
                   stl::RetrieveKey());
    graph::indexedGraph putative_graph(set_viewIds, resulting_pairs);
    graph::exportToGraphvizData(
        stlplus::create_filespec(base_folder, "retrieval_pairs"), putative_graph);
  }
}

#include <Eigen/Dense>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>

namespace cereal {
template <class Archive, class Derived>
inline typename std::enable_if<
    traits::is_output_serializable<BinaryData<typename Derived::Scalar>,
                                   Archive>::value, void>::type
save(Archive &ar, Eigen::PlainObjectBase<Derived> const &m) {
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  if (ArrT::RowsAtCompileTime == Eigen::Dynamic) ar(m.rows());
  if (ArrT::ColsAtCompileTime == Eigen::Dynamic) ar(m.cols());
  ar(binary_data(m.data(), m.size() * sizeof(typename Derived::Scalar)));
}

template <class Archive, class Derived>
inline typename std::enable_if<
    traits::is_input_serializable<BinaryData<typename Derived::Scalar>,
                                  Archive>::value, void>::type
load(Archive &ar, Eigen::PlainObjectBase<Derived> &m) {
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  Eigen::Index rows = ArrT::RowsAtCompileTime, cols = ArrT::ColsAtCompileTime;
  if (rows == Eigen::Dynamic) ar(rows);
  if (cols == Eigen::Dynamic) ar(cols);
  m.resize(rows, cols);
  ar(binary_data(m.data(),
                 static_cast<std::size_t>(rows * cols *
                                          sizeof(typename Derived::Scalar))));
}
}  // namespace cereal

int main(int argc, char **argv) {
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDirectory = "";
  std::string sPairFile = "vlad_pairs.txt";
  int32_t num_neighbors = 0;
  int32_t codebook_size = 128;
  int32_t vlad_flavor =
      static_cast<int>(VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW);
  int32_t max_feats = -1;
  uint32_t ui_max_cache_size = 0;

  // required
  cmd.add(make_option('i', sSfM_Data_Filename, "input_file"));
  cmd.add(make_option('o', sMatchesDirectory, "out_dir"));

  // optional
  cmd.add(make_option('n', num_neighbors, "num_neighbors"));
  cmd.add(make_option('d', codebook_size, "codebook_size"));
  cmd.add(make_option('p', sPairFile, "pair_file"));
  cmd.add(make_option('v', vlad_flavor, "vlad_flavor"));
  cmd.add(make_option('c', ui_max_cache_size, "cache_size"));
  cmd.add(make_option('m', max_feats, "max_feats"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string &s) {
    std::cerr
        << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] a SfM_Data file\n"
        << "[-o|--out_dir path] output path where computed features are "
           "stored\n"
        << "[-p|--pair_file] name of the output pair file (def. "
           "vlad_pairs.txt)\n"
        << "[-n|--num_neighbors] num neighbors per image (<= 0: auto i.e. 30% "
           "of the whole set)\n"
        << "[-d|--codebook_size] size of the codebook (number of kmeans "
           "centroids) used to compute descriptor (default=128)\n"
        << "[-v|--vlad_flavor] VLAD flavor (default=" << vlad_flavor << "):\n"
        << "\t" << static_cast<int>(VLAD_NORMALIZATION::SIGNED_SQUARE_ROOTING)
        << ": \"Aggregating local descriptors into compact codes\". H. Jegou "
           "et al. PAMI 2012\n"
        << "\t" << static_cast<int>(VLAD_NORMALIZATION::INTRA_NORMALIZATION)
        << ": \"All About VLAD\". R. Arandjelovic and A. Zisserman. CVPR 2013\n"
        << "\t"
        << static_cast<int>(VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW)
        << ": \"Revisiting the VLAD image representation\". J. Delhumeau et "
           "al. ACM Multimedia 2013. \n"
        << "[-m|--max_feats] Max number of features to perform the learning "
           "step (<= 0: whole feature set is used, default="
        << max_feats << ")\n"
        << "[-c|--cache_size] Use a regions cache (only cache_size regions "
           "will be stored in memory)\n"
        << "\t"
        << "If not used, all regions will be loaded in memory." << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << " You called : "
            << "\n"
            << argv[0] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--out_dir " << sMatchesDirectory << "\n"
            << "--pair_file " << sPairFile << "\n"
            << "--num_neighbors " << num_neighbors << "\n"
            << "--codebook_size " << codebook_size << "\n"
            << "--vlad_flavor " << vlad_flavor << "\n"
            << "--max_feats " << max_feats << "\n"
            << std::endl;

  if (sMatchesDirectory.empty() || !stlplus::is_folder(sMatchesDirectory)) {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS))) {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfM_Data_Filename
              << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions
  // extraction)
  using namespace openMVG::features;
  const std::string sImage_describer =
      stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type =
      Init_region_type_from_file(sImage_describer);
  if (!regions_type) {
    std::cerr << "Invalid: " << sImage_describer << " regions type file."
              << std::endl;
    return EXIT_FAILURE;
  }

  // Load the corresponding view regions - for learning -
  std::shared_ptr<Regions_Provider> learning_regions_provider;

  if (max_feats <= 0) {
    if (ui_max_cache_size == 0) {
      // Default regions provider (load & store all regions in memory)
      learning_regions_provider = std::make_shared<Regions_Provider>();
    } else {
      // cached region provider (progressive loading)
      learning_regions_provider =
          std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
    }
  } else {
    learning_regions_provider = std::make_shared<Preemptive_Regions_Provider>(
        max_feats / sfm_data.GetViews().size());
  }

  system::LoggerProgress progress;
  if (!learning_regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &progress))
  {
    std::cerr << std::endl
              << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  // Default parameters for num_neighbors
  if (num_neighbors <= 0) {
    num_neighbors = static_cast<int>(std::ceil(sfm_data.views.size() * 0.3));
  }

  if (num_neighbors > sfm_data.views.size()) {
    num_neighbors = sfm_data.views.size() - 1;
  }

  const VLAD_NORMALIZATION vlad_normalization =
      static_cast<VLAD_NORMALIZATION>(vlad_flavor);

  const size_t base_descriptor_length =
      learning_regions_provider->getRegionsType()->DescriptorLength();
  const size_t vlad_descriptor_length = base_descriptor_length * codebook_size;

  // -----------------------------
  // VLAD Computation
  // Learn the codebook: compute K centroids
  // Embedding: compute descriptors to centroid ids
  // Aggregation: Sum descriptor-to-centroid distances
  // Normalization: Each VLAD flavor come with different schemes of
  // normalization (mainly to handle visual bursts)
  // -----------------------------
  
  std::unique_ptr<VLADBase> vlad_builder;
  if (dynamic_cast<const SIFT_Regions*>(regions_type.get())){
    OPENMVG_LOG_INFO << "SIFT";
    vlad_builder.reset(new VLAD<SIFT_Regions>); 
  }
  else
  if (dynamic_cast<const AKAZE_Float_Regions*>(regions_type.get())) {
    OPENMVG_LOG_INFO << "AKAZE";
    vlad_builder.reset(new VLAD<AKAZE_Float_Regions>);
  }
  else {
    OPENMVG_LOG_ERROR << "VLAD does not support this Regions type.";
    OPENMVG_LOG_ERROR << "Please consider add the specialization here.";
    return EXIT_FAILURE;
  }

  std::vector<IndexT> view_ids;
  for (const auto &view : sfm_data.GetViews()) {
    const auto &view_id = view.first;
    view_ids.push_back(view_id);
  }

  // Convert input regions to array
  VLADBase::DescriptorVector descriptor_array = vlad_builder->RegionsToCodebook(
    view_ids,
    learning_regions_provider);

  std::cout << "Using # features for learning: " << descriptor_array.size()
          << std::endl;

  VLADBase::DescriptorVector codebook;
  codebook = vlad_builder->BuildCodebook(descriptor_array);
  
  // Freeing some memory
  descriptor_array.clear();
  descriptor_array.shrink_to_fit();

  std::unique_ptr<features::Regions> codebook_regions(regions_type->EmptyClone());
  vlad_builder->CodebookToRegions(codebook_regions, codebook);
  
  std::shared_ptr<Regions_Provider> embedding_regions_provider;
  if (max_feats <= 0) {
    embedding_regions_provider = std::move(learning_regions_provider);
  } else {
    // cached region provider (progressive loading)
    if (ui_max_cache_size == 0) {
      embedding_regions_provider = std::make_shared<Regions_Provider>();
    } else {
      embedding_regions_provider =
          std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
    }
    if (!embedding_regions_provider->load(sfm_data, sMatchesDirectory,
                                          regions_type, &progress)) {
      std::cerr << std::endl << "Invalid regions." << std::endl;
      return EXIT_FAILURE;
    }
  }

  VLADBase::VladMatrixType vlad_image_descriptors =
    vlad_builder->ComputeVLADEmbedding(
      view_ids,
      codebook_regions,
      embedding_regions_provider);

  // release the region provider
  embedding_regions_provider.reset();

  // Data structures to store the Results
  Pair_Set resulting_pairs;
  using DescendingIndexedPairwiseSimilarity =
      IndexedPairwiseSimilarity<std::greater<double>>;
  DescendingIndexedPairwiseSimilarity result_ordered_by_similarity;

  //
  // Retrieval
  //
  progress.Restart(
      sfm_data.GetViews().size(), "- VLAD Retrieval... -");
  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
  #endif
  for (int i = 0; i < sfm_data.GetViews().size(); i++) {
    const auto &view_id = sfm_data.GetViews().at(i)->id_view;
    const VLADBase::VladInternalType *query =
        vlad_image_descriptors.col(view_id).data();

    matching::ArrayMatcherBruteForce<VLADBase::VladInternalType,
                                    matching::LInner<VLADBase::VladInternalType>>
        matcher;

    if (!matcher.Build(vlad_image_descriptors.data(), sfm_data.views.size(),
                      vlad_descriptor_length)) {
      std::cout << "Error::Build" << std::endl;
    }

    const size_t NN = num_neighbors + 1;  // num_neighbors + 1 (the query vector
                                          // itself is part of the database)
    IndMatches nearest_neighbor_ids;
    std::vector<VLADBase::VladInternalType> nearest_neighbor_similarities;
    if (!matcher.SearchNeighbours(query, 1, &nearest_neighbor_ids,
                                  &nearest_neighbor_similarities, NN))
      std::cout << "Error::SearchNeighbours" << std::endl;
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
    #endif
    {
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
    #endif
    {
      for (int id = 0; id < nearest_neighbor_ids.size(); ++id) {
        const auto pair = nearest_neighbor_ids[id];
        if (view_id == pair.j_) continue;  // Ignore if we find the same image
        const auto similarity = -1. * nearest_neighbor_similarities[id];
        resulting_pairs.insert(
            {std::min(view_id, pair.j_), std::max(view_id, pair.j_)});
        result_ordered_by_similarity[view_id].insert({similarity, pair.j_});
      }
      ++progress;
    }
    }
  }

  // Export pairs into a text file
  savePairs(sPairFile, resulting_pairs);

  saveAdjacencyMatrixViewGraph(resulting_pairs, sfm_data, stlplus::basename_part(sMatchesDirectory));

  // Export the retrieval matrix
  saveRetrievalMatrix(
      stlplus::create_filespec(sMatchesDirectory, "retrieval_matches_matrix.svg"),
      sfm_data, result_ordered_by_similarity);

  // Export the sim file
  std::string sSimFile = stlplus::create_filespec(
      stlplus::folder_part(sPairFile),
      stlplus::create_filename(stlplus::basename_part(sPairFile) + "_simscores",
                              ".txt"));

  savePairwiseSimilarityScores(sSimFile, result_ordered_by_similarity);

  return EXIT_SUCCESS;  
}
