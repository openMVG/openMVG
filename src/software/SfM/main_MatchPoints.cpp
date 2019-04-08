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
               vector<int> &good_points, vector<int> &views_indices)
{
  // load sfm_data
  if (!Load(sfm_data, sfm_data_file, ESfM_Data(VIEWS|STRUCTURE)))
  {
    cerr << "Failed to load " << sfm_data_file << endl;
    return false;
  }

  // load descs, feats and featmap from file
  map<int, DescsT> view_all_descs;
  map<int, std::vector<Rich_SIOPointFeature>> view_all_feats;
  map<int, map<int/*old_feat*/, int/*new_feat*/>> view_featmap;
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
          view_all_descs[iter->first] = std::move(descs);
          view_all_feats[iter->first] = std::move(feats);
          view_featmap[iter->first] = std::move(featmap);
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
      if (view_featmap.find(ob.first) != view_featmap.end())
      {
        const auto &featmap = view_featmap.at(ob.first);
        if (featmap.find(ob.second.id_feat) != featmap.end())
        {
          const auto &descs = view_all_descs.at(ob.first);
          const auto &feats = view_all_feats.at(ob.first);
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
  if (!load_data(argv[1], argv[3], output_dir, "front", front_sfm_data, front_descs, front_feats, front_points, front_views))
    return 1;
  // back
  SfM_Data back_sfm_data;
  DescsT back_descs;
  vector<Rich_SIOPointFeature> back_feats;
  vector<int> back_points;
  vector<int> back_views;
  if (!load_data(argv[2], argv[3], output_dir, "back", back_sfm_data, back_descs, back_feats, back_points, back_views))
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

  IndMatches vec_putatives_matches;
  vec_putatives_matches.reserve(vec_Indice.size());
  MetricT metric;
  for (const auto &pair : vec_Indice)
  {
    if (metric(front_descs[pair.i_].data(), back_descs[pair.j_].data(), DescriptorT::static_size) < 2 * 255 * 255)
    {
        vec_putatives_matches.emplace_back(pair);
    }
  }

  // for each front point, count its descs matching result
  map<int, vector<int>> front_point_match_result;
#ifdef OUTPUT_DEBUG_FILE
  map<IntPair, vector<float>> feat_match_data;
#endif
  if (!vec_putatives_matches.empty())
  {
    for (const auto &putative_match : vec_putatives_matches)
    {
      int fPIndex = front_points[putative_match.i_];
      int bPIndex = back_points[putative_match.j_];
      front_point_match_result[fPIndex].emplace_back(bPIndex);

#ifdef OUTPUT_DEBUG_FILE
      auto f_feat = front_feats[putative_match.i_];
      auto b_feat = back_feats[putative_match.j_];
      auto &data_vec = feat_match_data[IntPair(front_views[putative_match.i_], back_views[putative_match.j_])];
      data_vec.emplace_back(f_feat.x());
      data_vec.emplace_back(f_feat.y());
      data_vec.emplace_back(b_feat.x());
      data_vec.emplace_back(b_feat.y());
#endif
    }
#ifdef OUTPUT_DEBUG_FILE
    using fm_pair = pair<IntPair, vector<float>>;
    vector<fm_pair> sort_feat_match_data(feat_match_data.begin(), feat_match_data.end());
    sort(sort_feat_match_data.begin(), sort_feat_match_data.end(), [](const fm_pair &a, const fm_pair &b) -> bool {
        return a.second.size() > b.second.size();
    });

    string feat_match = create_filespec(output_dir, "matched_feat_result", "txt");
    ofstream feat_match_point(feat_match, std::ios::out);
    for (const auto &match_data : sort_feat_match_data)
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
  string outputPointsMatch = create_filespec(output_dir, "matched_point_result", "txt");
  ofstream points_match(outputPointsMatch, std::ios::out);
  for (auto& m : front_point_match_result)
  {
    // find the best match for each front point
    int front_point = m.first;
    vector<int> &data_vec = m.second;
    sort(data_vec.begin(), data_vec.end());

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

    // output the besh match result: front index, back index, score, front point's pos, back point's pos, score / descs count
    const auto &frontP = front_sfm_data.GetLandmarks().at(front_point);
    const auto &backP = back_sfm_data.GetLandmarks().at(best_match_back_point);
    points_match << front_point << " " << best_match_back_point << " " << best_match_score << " " << data_vec.size() << endl;
    points_match << frontP.X[0] << " " << frontP.X[1] << " " << frontP.X[2] << " "
                 << backP.X[0] << " " << backP.X[1] << " " << backP.X[2] << endl;

    // output feats pos in each image
    for (const auto &ob : frontP.obs)
    {
      points_match << front_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                   << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match << endl;

    // output feats pos in each image
    for (const auto &ob : backP.obs)
    {
      points_match << back_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                   << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match << endl;
  }
  points_match.close();

  return 0;
}
