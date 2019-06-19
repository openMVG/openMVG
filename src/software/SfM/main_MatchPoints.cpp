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

#define OUTPUT_DEBUG_FILE

using DescriptorT = Descriptor<unsigned char, 128>;
using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;
using IntPair = pair<int, int>;

bool load_featmap(const string &featmap_file, map<int, int> &featmap)
{
  ifstream stream(featmap_file, std::ios::in);
  if (!stream.is_open())
  {
      return false;
  }
  int old_feat_index = 0;
  int new_feat_index = 0;
  while (stream >> old_feat_index)
  {
    featmap[old_feat_index] = new_feat_index++;
  }
  stream.close();
  return true;
}

bool load_data(const string &sfm_data_file, const string &matches_dir, const string &output_dir, const string &prefix,
               SfM_Data &sfm_data, DescsT &good_descs, vector<Rich_SIOPointFeature> &good_feats,
               vector<int> &good_points, vector<int> &views_indices, map<int, map<int, int>> &view_id2featmap)
{
  // load sfm_data
  if (!Load(sfm_data, sfm_data_file, ESfM_Data(VIEWS|STRUCTURE)))
  {
    cerr << "Failed to load " << sfm_data_file << endl;
    return false;
  }

  good_descs.clear();
  good_feats.clear();
  good_points.clear();
  views_indices.clear();
  view_id2featmap.clear();

  // load descs, feats and featmap from file
  map<int, DescsT> view_id2descs;
  map<int, std::vector<Rich_SIOPointFeature>> view_id2feats;
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
  for (auto iter = sfm_data.GetViews().begin(); iter != sfm_data.GetViews().end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
    {
      const string sImageName = iter->second->s_Img_path;
      const string basename = basename_part(sImageName);
      const string descFile = create_filespec(matches_dir, basename, ".desc");
      const string featFile = create_filespec(matches_dir, basename, "feat");
      const string featmapFile = create_filespec(matches_dir, basename, "featmap");

      DescsT descs;
      std::vector<Rich_SIOPointFeature> feats;
      map<int, int> featmap;
      if (loadDescsFromBinFile(descFile, descs) &&
          loadFeatsFromFile(featFile, feats) &&
          load_featmap(featmapFile, featmap))
      {
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          view_id2descs[iter->first] = std::move(descs);
          view_id2feats[iter->first] = std::move(feats);
          view_id2featmap[iter->first] = std::move(featmap);
        }
      }
    }
  }

  // keep data which relates to 3d points
#ifdef OUTPUT_DEBUG_FILE
  string outputFilteredFeats = create_filespec(output_dir, string("filtered_feat_") + prefix, "txt");
  ofstream stream_feats(outputFilteredFeats, std::ios::out);
  string outputFilteredPoints = create_filespec(output_dir, string("filtered_point_") + prefix, "txt");
  ofstream stream_points(outputFilteredPoints, std::ios::out);
  set<int> points_set;
#endif
  for (const auto &structure : sfm_data.GetLandmarks())
  {
    for (const auto &ob : structure.second.obs)
    {
      if (view_id2featmap.find(ob.first) != view_id2featmap.end())
      {
        const auto &featmap = view_id2featmap.at(ob.first);
        if (featmap.find(ob.second.id_feat) != featmap.end())
        {
          const auto &descs = view_id2descs.at(ob.first);
          const auto &feats = view_id2feats.at(ob.first);
          int new_feat_index = featmap.at(ob.second.id_feat);
          good_descs.emplace_back(descs[new_feat_index]);
          good_feats.emplace_back(feats[new_feat_index]);
          good_points.emplace_back(structure.first);
          views_indices.emplace_back(ob.first);
#ifdef OUTPUT_DEBUG_FILE
          stream_feats << sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                       << feats[new_feat_index].x() << " " << feats[new_feat_index].y() << endl;
          if (points_set.count(structure.first) == 0)
          {
              stream_points << structure.first << " " << structure.second.X[0] << " "
                            << structure.second.X[1] << " " << structure.second.X[2] << endl;
          }
          points_set.emplace(structure.first);
#endif
        }
      }
    }
  }
#ifdef OUTPUT_DEBUG_FILE
  stream_feats.close();
  stream_points.close();
#endif

  return true;
}

int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    cerr << "Usage: openMVG_main_MatchPoints front_sfm_data_all back_sfm_data_all merge_matches_dir output_dir" << endl;
    return 1;
  }
  string output_dir = argv[4];

  // load sfm_data_all, build good feat info: descs, feats, point map, view map
  // front
  SfM_Data front_sfm_data;
  DescsT front_descs;
  vector<Rich_SIOPointFeature> front_feats;
  vector<int> front_points;
  vector<int> front_views;
  map<int, map<int/*old_feat*/, int/*new_feat*/>> front_view_id2featmap;
  if (!load_data(argv[1], argv[3], output_dir, "front", front_sfm_data, front_descs, front_feats, front_points, front_views, front_view_id2featmap))
    return 1;
  // back
  SfM_Data back_sfm_data;
  DescsT back_descs;
  vector<Rich_SIOPointFeature> back_feats;
  vector<int> back_points;
  vector<int> back_views;
  map<int, map<int/*old_feat*/, int/*new_feat*/>> back_view_id2featmap;
  if (!load_data(argv[2], argv[3], output_dir, "back", back_sfm_data, back_descs, back_feats, back_points, back_views, back_view_id2featmap))
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
  IndMatches vec_Match;
  std::vector<DistanceType> vec_Distance;
  // Search NNN__ neighbours for each query descriptor
  matcher.SearchNeighbours(queries, front_descs.size(), &vec_Match, &vec_Distance, NNN__);

  IndMatches vec_putative_match;
  vec_putative_match.reserve(vec_Match.size());
  MetricT metric;
  for (const auto &match : vec_Match)
  {
    if (metric(front_descs[match.i_].data(), back_descs[match.j_].data(), DescriptorT::static_size) < 2 * 255 * 255)
    {
        vec_putative_match.emplace_back(match);
    }
  }

  // for each front point, record its descs matching result
  map<int, vector<int>> point_match_info;
#ifdef OUTPUT_DEBUG_FILE
  map<IntPair, vector<float>> feat_match_info;
#endif
  if (!vec_putative_match.empty())
  {
    for (const auto &putative_match : vec_putative_match)
    {
      int fp_index = front_points[putative_match.i_];
      int bp_index = back_points[putative_match.j_];
      point_match_info[fp_index].emplace_back(bp_index);

#ifdef OUTPUT_DEBUG_FILE
      auto fv_index = front_views[putative_match.i_];
      auto bv_index = back_views[putative_match.j_];
      auto &match_vec = feat_match_info[IntPair(fv_index, bv_index)];

      auto f_feat = front_feats[putative_match.i_];
      auto b_feat = back_feats[putative_match.j_];
      match_vec.emplace_back(f_feat.x());
      match_vec.emplace_back(f_feat.y());
      match_vec.emplace_back(b_feat.x());
      match_vec.emplace_back(b_feat.y());
#endif
    }
#ifdef OUTPUT_DEBUG_FILE
    using fm_pair = pair<IntPair, vector<float>>;
    vector<fm_pair> sorted_feat_match_info(feat_match_info.begin(), feat_match_info.end());
    sort(sorted_feat_match_info.begin(), sorted_feat_match_info.end(), [](const fm_pair &a, const fm_pair &b) -> bool {
        return a.second.size() > b.second.size();
    });

    string feat_match_path = create_filespec(output_dir, "matched_feat_result", "txt");
    ofstream feat_match_result(feat_match_path, std::ios::out);
    for (const auto &match_info : sorted_feat_match_info)
    {
      feat_match_result << front_sfm_data.GetViews().at(match_info.first.first).get()->s_Img_path << " "
                        << back_sfm_data.GetViews().at(match_info.first.second).get()->s_Img_path << endl;
      for (const auto &feat_pos : match_info.second)
      {
        feat_match_result << feat_pos << " ";
      }
      feat_match_result << endl;
    }
    feat_match_result.close();
#endif
  }

  if (point_match_info.size() == 0) {
    cerr << "Failed to find any point match" << endl;
    return 1;
  }

  using pm_info = tuple<int, int, int, int>; // fp index, bp index, score, valid descs count
  vector<pm_info> sorted_point_match_info;
  sorted_point_match_info.reserve(point_match_info.size());
  for (auto& m : point_match_info)
  {
    // for each front point, find its best match
    int fp_index = m.first;
    vector<int> &match_vec = m.second;
    // find the most frequent element in a vector
    sort(match_vec.begin(), match_vec.end());
    int best_score = 0;
    int best_bp_index = -1;
    for (int j = 0; j < match_vec.size(); ++j)
    {
      int bp_index = match_vec[j];
      int count = 1;
      for (int k = j + 1; k < match_vec.size() && match_vec[k] == bp_index; ++k)
      {
        ++count;
        ++j;
      }
      if (count > best_score)
      {
        best_score = count;
        best_bp_index = bp_index;
      }
    }

    sorted_point_match_info.emplace_back(fp_index, best_bp_index, best_score, match_vec.size());
  }
  sort(sorted_point_match_info.begin(), sorted_point_match_info.end(), [](const pm_info &a, const pm_info &b) -> bool {
        return std::get<2>(a) > std::get<2>(b);
    });

  // output the best match for each front point
  string point_match_file = create_filespec(output_dir, "matched_point_result", "txt");
  ofstream points_match_result(point_match_file, std::ios::out);
  for (auto& m : sorted_point_match_info)
  {
    auto fp_index = std::get<0>(m);
    auto bp_index = std::get<1>(m);
    auto score = std::get<2>(m);
    auto vote_num = std::get<3>(m);

    const auto &frontP = front_sfm_data.GetLandmarks().at(fp_index);
    const auto &backP = back_sfm_data.GetLandmarks().at(bp_index);
    int back_good_feat_num = 0;
    for (const auto &ob : backP.obs)
    {
      const auto &featmap = back_view_id2featmap[ob.first];
      if (featmap.find(ob.second.id_feat) != featmap.end())
        ++back_good_feat_num;
    }
    // if 3d point do not have enough features after modelmask filter
    // its 3d position is probabily not that good
    if (vote_num < 2 || back_good_feat_num < 2)
      continue;

    // output fp index, bp index, score, valid descs count
    points_match_result << fp_index << " " << bp_index << " " << score << " " << vote_num << endl;
    // output fp pos, bp pos
    points_match_result << frontP.X[0] << " " << frontP.X[1] << " " << frontP.X[2] << " "
                        << backP.X[0] << " " << backP.X[1] << " " << backP.X[2] << endl;

    // output feats pos in each view
    for (const auto &ob : frontP.obs)
    {
      const auto &featmap = front_view_id2featmap[ob.first];
      if (featmap.find(ob.second.id_feat) != featmap.end())
        points_match_result << front_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                            << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match_result << endl;
    for (const auto &ob : backP.obs)
    {
      const auto &featmap = back_view_id2featmap[ob.first];
      if (featmap.find(ob.second.id_feat) != featmap.end())
        points_match_result << back_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                            << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match_result << endl;
  }
  points_match_result.close();

  return 0;
}
