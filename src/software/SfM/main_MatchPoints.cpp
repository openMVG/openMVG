#include <cereal/archives/json.hpp>

#include "nonFree/rich_sift/Rich_SIFT_describer_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/numeric/numeric.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace std;
using namespace stlplus;
using namespace openMVG::sfm;
using namespace openMVG::features;
using namespace openMVG::matching;

#define OUTPUT_TEMP_FILE
#ifdef OUTPUT_TEMP_FILE
const string OUTPUT_TEST_DIR = "/home/dushuai/work/test/match_point_result";
#endif

using DescriptorT = Descriptor<unsigned char, 128>;
using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;
using IntPair = pair<int, int>;

bool load_featmap(const string &featmap_file, map<int, int> &featmap)
{
  ifstream stream(featmap_file, std::ios::in);
  int new_feat_index = 0;
  while (!stream.eof())
  {
    int old_feat_index;
    stream >> old_feat_index;
    featmap[old_feat_index] = new_feat_index++;
  }
  stream.close();
}

bool load_data(const string &data_dir, const string &prefix,
               SfM_Data &sfm_data, DescsT &good_descs, vector<Rich_SIOPointFeature> &good_feats,
               vector<int> &points, vector<int> &front_views)
{
  string sfm_data_file = create_filespec(data_dir, prefix + "/output/reconstruction_sequential/sfm_data_all", "json");
  string matches_dir = create_filespec(data_dir, "merge/output/matches");

  // load sfm_data
  if (!Load(sfm_data, sfm_data_file, ESfM_Data(VIEWS|STRUCTURE)))
  {
    cerr << "Failed to load " << sfm_data_file << endl;
    return false;
  }

  // load descs and featmap from file
  map<int, DescsT> view_all_descs;
  map<int, map<int/*old_feat*/, int/*new_feat*/>> view_featmap;
  map<int, std::vector<Rich_SIOPointFeature>> view_all_feats;

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
  for (auto iter = sfm_data.GetViews().begin(); iter != sfm_data.GetViews().end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
    {
      const string sImageName = create_filespec(sfm_data.s_root_path, iter->second->s_Img_path);
      const string basename = basename_part(sImageName);
      const string descFile = create_filespec(matches_dir, basename, ".desc");
      const string featFile = create_filespec(matches_dir, basename, "feat");
      const string featmapFile = create_filespec(matches_dir, basename, "featmap");

      DescsT descs;
      map<int, int> featmap;
      std::vector<Rich_SIOPointFeature> feats;
      if (loadDescsFromBinFile(descFile, descs) &&
          load_featmap(featmapFile, featmap) &&
          loadFeatsFromFile(featFile, feats))
      {
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          view_all_descs[iter->first] = std::move(descs);
          view_featmap[iter->first] = std::move(featmap);
          view_all_feats[iter->first] = std::move(feats);
        }
      }
    }
  }

  // build view_good_descs and view_points
#ifdef OUTPUT_TEMP_FILE
  string outputFilteredFeats = create_filespec(OUTPUT_TEST_DIR, prefix + "_filtered_feat", "txt");
  ofstream stream(outputFilteredFeats, std::ios::out);
#endif
  for (const auto &structure : sfm_data.GetLandmarks())
  {
    for (const auto &ob : structure.second.obs)
    {
      if (view_featmap.find(ob.first) != view_featmap.end())
      {
        const auto &featmap = view_featmap.at(ob.first);
        if (featmap.find(ob.second.id_feat) != featmap.end())
        {
          int new_feat_index = featmap.at(ob.second.id_feat);
          const auto &descs = view_all_descs.at(ob.first);
          const auto &feats = view_all_feats.at(ob.first);
          good_descs.emplace_back(descs[new_feat_index]);
          good_feats.emplace_back(feats[new_feat_index]);
          points.emplace_back(structure.first);
          front_views.emplace_back(ob.first);
#ifdef OUTPUT_TEMP_FILE
          stream << sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                 << feats[new_feat_index].x() << " " << feats[new_feat_index].y() << endl;
#endif
        }
      }
    }
  }
#ifdef OUTPUT_TEMP_FILE
  stream.close();
#endif

  return true;
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    cerr << "Usage: openMVG_main_MatchPoints root_data_dir output_dir" << endl;
    return 1;
  }

  string root_data_dir = argv[1];

  // load sfm_data_all, build good view_descs and view_points
  // front
  SfM_Data front_sfm_data;
  DescsT front_descs;
  vector<Rich_SIOPointFeature> front_feats;
  vector<int> front_points;
  vector<int> front_views;
  if (!load_data(root_data_dir, "front", front_sfm_data, front_descs, front_feats, front_points, front_views))
    return 1;
  // back
  SfM_Data back_sfm_data;
  DescsT back_descs;
  vector<Rich_SIOPointFeature> back_feats;
  vector<int> back_points;
  vector<int> back_views;
  if (!load_data(root_data_dir, "back", back_sfm_data, back_descs, back_feats, back_points, back_views))
    return 1;


  // match descs
  using MetricT = L2<unsigned char>;
  using ArrayMatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
  using Scalar = typename ArrayMatcherT::ScalarT;
  using DistanceType = typename ArrayMatcherT::DistanceType;

  ArrayMatcherT matcher;
  matcher.Build(reinterpret_cast<const Scalar*>(&back_descs[0]), back_descs.size(), 128);

  auto queries = reinterpret_cast<const Scalar*>(&front_descs[0]);
  const size_t NNN__ = 1;
  IndMatches vec_Indice;
  std::vector<DistanceType> vec_Distance;
  // Search NNN__ neighbours for each query descriptor
  matcher.SearchNeighbours(queries, front_descs.size(), &vec_Indice, &vec_Distance, NNN__);
  IndMatches vec_putatives_matches = vec_Indice;

  map<int, vector<int>> front_point_match_result;
#ifdef OUTPUT_TEMP_FILE
  map<IntPair, vector<float>> feat_match_data;
#endif
  if (!vec_putatives_matches.empty())
  {
    for (const auto &putative_match : vec_putatives_matches)
    {
      int fPIndex = front_points[putative_match.i_];
      int bPIndex = back_points[putative_match.j_];
      front_point_match_result[fPIndex].emplace_back(bPIndex);

#ifdef OUTPUT_TEMP_FILE
      auto f_feat = front_feats[putative_match.i_];
      auto b_feat = back_feats[putative_match.j_];
      auto &data_vec = feat_match_data[IntPair(front_views[putative_match.i_], back_views[putative_match.j_])];
      data_vec.emplace_back(f_feat.x());
      data_vec.emplace_back(f_feat.y());
      data_vec.emplace_back(b_feat.x());
      data_vec.emplace_back(b_feat.y());
#endif
    }
#ifdef OUTPUT_TEMP_FILE
    string feat_match = create_filespec(OUTPUT_TEST_DIR, "match_feat_3d_points_result", "txt");
    ofstream feat_match_point(feat_match, std::ios::out);
    for (const auto &match_data : feat_match_data)
    {
      feat_match_point << front_sfm_data.GetViews().at(match_data.first.first).get()->s_Img_path << " "
                       << back_sfm_data.GetViews().at(match_data.first.second).get()->s_Img_path << endl;
      for (const auto &pos : match_data.second)
      {
        feat_match_point << pos << " ";
      }
      feat_match_point << endl;
    }
    feat_match_point.close();
#endif
  }

  // output the best match for each front point
  string outputPointsMatch = create_filespec(argv[2], "match_point_result", "txt");
  ofstream points_match(outputPointsMatch, std::ios::out);
  for (auto& m : front_point_match_result)
  {
    int front_point = m.first;
    vector<int> &data_vec = m.second;
    sort(data_vec.begin(), data_vec.end());

    // find the best match
    int best_match_score = 0;
    int best_match_back_point = -1;
    for (int j = 0; j < data_vec.size(); ++j)
    {
      int score = 1;
      int index = data_vec[j];
      for (int k = j + 1; k < data_vec.size() && data_vec[k] == index; ++k)
      {
        ++score;
        ++j;
      }
      if (score > best_match_score)
      {
        best_match_score = score;
        best_match_back_point = index;
      }
    }

    const auto &frontP = front_sfm_data.GetLandmarks().at(front_point);
    const auto &backP = back_sfm_data.GetLandmarks().at(best_match_back_point);
    points_match << front_point << " " << best_match_back_point << " " << best_match_score << " "
                 << frontP.X[0] << " " << frontP.X[1] << " " << frontP.X[2] << " "
                 << backP.X[0] << " " << backP.X[1] << " " << backP.X[2] << " " << (float)best_match_score / data_vec.size() << endl;

    for (const auto &ob : frontP.obs)
    {
      points_match << front_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                   << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match << endl;

    for (const auto &ob : backP.obs)
    {
      points_match << back_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                   << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match << endl;
  }
  points_match.close();

  // remove redundant points
  string frontSparsePointsFile = create_filespec(argv[2], "front_point_cloud", "txt");
  ofstream front_point_cloud(frontSparsePointsFile, std::ios::out);
  set<int> front_good_points;
  front_good_points.insert(front_points.begin(), front_points.end());
  for (const auto &lm : front_sfm_data.GetLandmarks())
  {
    if (front_good_points.count(lm.first) != 0)
    {
      front_point_cloud << lm.first << " " << lm.second.X[0] << " " << lm.second.X[1] << " " << lm.second.X[2] << endl;
    }
  }
  front_point_cloud.close();

  string backSparsePointsFile = create_filespec(argv[2], "back_point_cloud", "txt");
  ofstream back_point_cloud(backSparsePointsFile, std::ios::out);
  set<int> back_good_points;
  back_good_points.insert(back_points.begin(), back_points.end());
  for (const auto &lm : back_sfm_data.GetLandmarks())
  {
    if (back_good_points.count(lm.first) != 0)
    {
      back_point_cloud << lm.first << " " << lm.second.X[0] << " " << lm.second.X[1] << " " << lm.second.X[2] << endl;
    }
  }
  back_point_cloud.close();

  return 1;
}






//#include <cereal/archives/json.hpp>
//
//#include "nonFree/rich_sift/Rich_SIFT_describer_io.hpp"
//#include "openMVG/sfm/sfm_data.hpp"
//#include "openMVG/sfm/sfm_data_io.hpp"
//#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
//#include "openMVG/features/feature.hpp"
//#include "openMVG/matching/matcher_cascade_hashing.hpp"
//#include "openMVG/matching/matcher_brute_force.hpp"
//#include "openMVG/matching/matching_filters.hpp"
//#include "openMVG/numeric/numeric.h"
//
//#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
//
//using namespace std;
//using namespace stlplus;
//using namespace openMVG::sfm;
//using namespace openMVG::features;
//using namespace openMVG::matching;
//
//#define OUTPUT_TEMP_FILE
//#ifdef OUTPUT_TEMP_FILE
//const string OUTPUT_TEST_DIR = "/home/dushuai/work/test/match_point_result";
//#endif
//
//using DescriptorT = Descriptor<unsigned char, 128>;
//using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;
//using IntPair = pair<int, int>;
//
//bool load_featmap(const string &featmap_file, map<int, int> &featmap)
//{
//  ifstream stream(featmap_file, std::ios::in);
//  int new_feat_index = 0;
//  while (!stream.eof())
//  {
//    int old_feat_index;
//    stream >> old_feat_index;
//    featmap[old_feat_index] = new_feat_index++;
//  }
//  stream.close();
//}
//
//bool load_data(const string &data_dir, const string &prefix,
//  SfM_Data &sfm_data, map<int, DescsT> &view_good_descs, map<int, vector<int>> &view_points)
//{
//  string sfm_data_file = create_filespec(data_dir, prefix + "/output/reconstruction_sequential/sfm_data_all", "json");
//  string matches_dir = create_filespec(data_dir, "merge/output/matches");
//
//  // load sfm_data
//  if (!Load(sfm_data, sfm_data_file, ESfM_Data(VIEWS|STRUCTURE)))
//  {
//    cerr << "Failed to load " << sfm_data_file << endl;
//    return false;
//  }
//
//  // load descs and featmap from file
//  map<int, DescsT> view_all_descs;
//  map<int, map<int/*old_feat*/, int/*new_feat*/>> view_featmap;
//#ifdef OUTPUT_TEMP_FILE
//  map<int, std::vector<Rich_SIOPointFeature>> view_all_feats;
//#endif
//
//#ifdef OPENMVG_USE_OPENMP
//#pragma omp parallel
//#endif
//  for (auto iter = sfm_data.GetViews().begin(); iter != sfm_data.GetViews().end(); ++iter)
//  {
//#ifdef OPENMVG_USE_OPENMP
//#pragma omp single nowait
//#endif
//    {
//      const string sImageName = create_filespec(sfm_data.s_root_path, iter->second->s_Img_path);
//      const string basename = basename_part(sImageName);
//      const string descFile = create_filespec(matches_dir, basename, ".desc");
//      const string featFile = create_filespec(matches_dir, basename, "feat");
//      const string featmapFile = create_filespec(matches_dir, basename, "featmap");
//
//      DescsT descs;
//      map<int, int> featmap;
//      std::vector<Rich_SIOPointFeature> feats;
//      if (loadDescsFromBinFile(descFile, descs) &&
//          load_featmap(featmapFile, featmap) &&
//          loadFeatsFromFile(featFile, feats))
//      {
//#ifdef OPENMVG_USE_OPENMP
//#pragma omp critical
//#endif
//        {
//          view_all_descs[iter->first] = std::move(descs);
//          view_featmap[iter->first] = std::move(featmap);
//#ifdef OUTPUT_TEMP_FILE
//          view_all_feats[iter->first] = std::move(feats);
//#endif
//        }
//      }
//    }
//  }
//
//  // build view_good_descs and view_points
//#ifdef OUTPUT_TEMP_FILE
//  string outputFilteredFeats = create_filespec(OUTPUT_TEST_DIR, prefix + "_filtered_feat", "txt");
//  ofstream stream(outputFilteredFeats, std::ios::out);
//#endif
//  for (const auto &structure : sfm_data.GetLandmarks())
//  {
//    for (const auto &ob : structure.second.obs)
//    {
//      if (view_featmap.find(ob.first) != view_featmap.end())
//      {
//        const auto &featmap = view_featmap.at(ob.first);
//        if (featmap.find(ob.second.id_feat) != featmap.end())
//        {
//          int new_feat_index = featmap.at(ob.second.id_feat);
//          const auto &descs = view_all_descs.at(ob.first);
//          view_good_descs[ob.first].emplace_back(descs[new_feat_index]);
//          view_points[ob.first].emplace_back(structure.first);
//#ifdef OUTPUT_TEMP_FILE
//          const auto &feats = view_all_feats.at(ob.first);
//          stream << sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
//                 << feats[new_feat_index].x() << " " << feats[new_feat_index].y() << endl;
//#endif
//        }
//      }
//    }
//  }
//#ifdef OUTPUT_TEMP_FILE
//  stream.close();
//#endif
//
//  return true;
//}
//
//int main(int argc, char *argv[])
//{
//  if (argc < 3)
//  {
//    cerr << "Usage: openMVG_main_MatchPoints root_data_dir output_dir" << endl;
//    return 1;
//  }
//
//  string root_data_dir = argv[1];
//
//  // load sfm_data_all, build good view_descs and view_points
//  // front
//  SfM_Data front_sfm_data;
//  map<int, DescsT> front_view_descs;
//  map<int, vector<int>> front_view_points;
//  if (!load_data(root_data_dir, "front", front_sfm_data, front_view_descs, front_view_points))
//    return 1;
//  // back
//  SfM_Data back_sfm_data;
//  map<int, DescsT> back_view_descs;
//  map<int, vector<int>> back_view_points;
//  if (!load_data(root_data_dir, "back", back_sfm_data, back_view_descs, back_view_points))
//    return 1;
//
//
//  // match descs
//  float fRatioDis = 0.8f;
//  map<IntPair, int> point_match_result;
//
//  using MetricT = L2<unsigned char>;
//  using ArrayMatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
//  using Scalar = typename ArrayMatcherT::ScalarT;
//  using DistanceType = typename ArrayMatcherT::DistanceType;
//
//#ifdef OUTPUT_TEMP_FILE
//  string feat_match = create_filespec(OUTPUT_TEST_DIR, "match_feat_3d_points_result", "txt");
//  ofstream feat_match_point(feat_match, std::ios::out);
//#endif
//  for (auto &back_ele : back_view_descs)
//  {
//    auto &back_descs = back_ele.second;
//
//    ArrayMatcherT matcher;
//    matcher.Build(reinterpret_cast<const Scalar*>(&back_descs[0]), back_descs.size(), 128);
//
//    for (const auto &front_ele : front_view_descs)
//    {
//      auto &front_descs = front_ele.second;
//
//      auto queries = reinterpret_cast<const Scalar*>(&front_descs[0]);
//      const size_t NNN__ = std::max<int>(back_descs.size() * 2 / 10, 2);
//      IndMatches vec_Indice;
//      std::vector<DistanceType> vec_Distance;
//      // Search NNN__ neighbours for each query descriptor
//      if (!matcher.SearchNeighbours(queries, front_descs.size(), &vec_Indice, &vec_Distance, NNN__))
//        continue;
//
//      // Filter the matches using a distance ratio test:
//      //   The probability that a match is correct is determined by taking
//      //   the ratio of distance from the closest neighbor to the distance
//      //   of the second closest.
//      std::vector<int> vec_nn_ratio_idx;
//      NNdistanceRatio(
//              vec_Distance.begin(), // distance start
//              vec_Distance.end(),   // distance end
//              NNN__, // Number of neighbor in iterator sequence (minimum required 2)
//              vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
//              openMVG::Square(fRatioDis));
//
//      IndMatches vec_putatives_matches;
//      vec_putatives_matches.reserve(vec_nn_ratio_idx.size());
//      for ( const auto & index : vec_nn_ratio_idx )
//      {
//        vec_putatives_matches.emplace_back(vec_Indice[index*NNN__]);
//      }
//      // Remove duplicates
//      IndMatch::getDeduplicated(vec_putatives_matches);
//
//      if (!vec_putatives_matches.empty())
//      {
//#ifdef OUTPUT_TEMP_FILE
//        feat_match_point << front_sfm_data.GetViews().at(front_ele.first).get()->s_Img_path << " "
//                         << back_sfm_data.GetViews().at(back_ele.first).get()->s_Img_path << endl;
//#endif
//        for (const auto &putative_match : vec_putatives_matches)
//        {
//          int fPIndex = front_view_points.at(front_ele.first)[putative_match.i_];
//          int bPIndex = back_view_points.at(back_ele.first)[putative_match.j_];
//          point_match_result[IntPair(fPIndex, bPIndex)]++;
//
//#ifdef OUTPUT_TEMP_FILE
//          auto f_feat_pos = front_sfm_data.GetLandmarks().at(fPIndex).obs.at(front_ele.first).x;
//          auto b_feat_pos = back_sfm_data.GetLandmarks().at(bPIndex).obs.at(back_ele.first).x;
//          feat_match_point << f_feat_pos.x() << " " << f_feat_pos.y() << " "
//                           << b_feat_pos.x() << " " << b_feat_pos.y() << " ";
//#endif
//        }
//#ifdef OUTPUT_TEMP_FILE
//        feat_match_point << endl;
//#endif
//      }
//    }
//  }
//#ifdef OUTPUT_TEMP_FILE
//  feat_match_point.close();
//#endif
//
//  // sort and output the top 100 matches for visualization
//  vector<pair<IntPair, int>> vec_matches;
//  vec_matches.reserve(point_match_result.size());
//  for (const auto &m : point_match_result)
//  {
//    vec_matches.emplace_back(pair<IntPair, int>(m.first, m.second));
//  }
//  sort(vec_matches.begin(), vec_matches.end(),
//       [](const pair<IntPair, int> &a, pair<IntPair, int> &b) -> bool {
//       return a.second > b.second;
//       });
//
//  string outputPointsMatch = create_filespec(argv[2], "match_point_result", "txt");
//  ofstream points_match(outputPointsMatch, std::ios::out);
//  for (int i = 0; i < point_match_result.size(); ++i)
//  {
//    const auto &pair = vec_matches[i].first;
//    const auto &frontP = front_sfm_data.GetLandmarks().at(pair.first);
//    const auto &backP = back_sfm_data.GetLandmarks().at(pair.second);
//    points_match << pair.first << " " << pair.second << " " << vec_matches[i].second << " "
//                 << frontP.X[0] << " " << frontP.X[1] << " " << frontP.X[2] << " "
//                 << backP.X[0] << " " << backP.X[1] << " " << backP.X[2] << endl;
//
//    for (const auto &ob : frontP.obs)
//    {
//      if (front_view_descs.find(ob.first) != front_view_descs.end())
//      {
//        points_match << front_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
//                     << ob.second.x[0] << " " << ob.second.x[1] << " ";
//      }
//    }
//    points_match << endl;
//
//    for (const auto &ob : backP.obs)
//    {
//      if (back_view_descs.find(ob.first) != back_view_descs.end())
//      {
//        points_match << back_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
//                     << ob.second.x[0] << " " << ob.second.x[1] << " ";
//      }
//    }
//    points_match << endl;
//  }
//  points_match.close();
//
//  // @TODO
//  // for each landmark in front, loop its features
//  //   for each feature, find the k-nearest features in matcher
//  // decide the best match for this landmark
//
//  return 1;
//}
