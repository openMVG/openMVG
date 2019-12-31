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

using DescriptorT = Descriptor<unsigned char, 128>;
using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;
using IntPair = pair<int, int>;

struct GoodFeatInfo {
  DescsT descs;
  vector<int> points;
  map<int, map<int/*old_feat*/, int/*new_feat*/>> featmap_map;
};

bool LoadFeatmap(const string &featmap_file, map<int, int> &featmap) {
  featmap.clear();
  ifstream stream(featmap_file, std::ios::in);
  if (!stream.is_open()) {
      return false;
  }
  int old_feat_index = 0;
  int new_feat_index = 0;
  while (stream >> old_feat_index) {
    featmap[old_feat_index] = new_feat_index++;
  }
  stream.close();
  return true;
}

bool LoadData(const SfM_Data &sfm_data, const string &matches_dir, const string &prefix,
  const string &output_dir, GoodFeatInfo &feat_info) {
  feat_info.descs.clear();
  feat_info.points.clear();
  feat_info.featmap_map.clear();

  // load descs and featmap from file
  map<int, DescsT> descs_map;
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel
#endif
  for (const auto &view : sfm_data.GetViews()) {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
    {
      string basename = basename_part(view.second->s_Img_path);

      DescsT descs;
      uint32_t descs_size = 0u;
      string descFile = create_filespec(matches_dir, basename, ".desc");
      if (loadDescsFromBinFile(descFile, descs)) {
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          descs_size = descs.size();
          descs_map[view.first] = std::move(descs);
        }
      }

      map<int, int> featmap;
      string featmapFile = create_filespec(matches_dir, basename, "featmap");
      if (LoadFeatmap(featmapFile, featmap)) {
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          feat_info.featmap_map[view.first] = std::move(featmap);
        }
      } else if (descs_size != 0) {
        for (uint32_t i = 0; i < descs_size; ++i) {
          featmap[i] = i;
        }
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          feat_info.featmap_map[view.first] = std::move(featmap);
        }
      }
    }
  }

  // keep data which relates to 3d points
  for (const auto &structure : sfm_data.GetLandmarks()) {
    for (const auto &ob : structure.second.obs) {
      auto iter1 = feat_info.featmap_map.find(ob.first);
      if (iter1 != feat_info.featmap_map.end()) {
        const auto &featmap = iter1->second;
        auto iter2 = featmap.find(ob.second.id_feat);
        if (iter2 != featmap.end()) {
          const auto &descs = descs_map.at(ob.first);
          int new_feat_index = iter2->second;
          feat_info.descs.emplace_back(descs[new_feat_index]);
          feat_info.points.emplace_back(structure.first);
        }
      }
    }
  }

  return true;
}

int main(int argc, char *argv[]) {
  if (argc < 5) {
    cerr << "Usage: openMVG_main_MatchPoints front_sfm_data_all back_sfm_data_all merge_matches_dir output_dir" << endl;
    return 1;
  }
  string front_sfm_data_file = argv[1];
  string back_sfm_data_file = argv[2];
  string output_dir = argv[4];

  // load sfm_data, build good feat info: descs, feats, point map, view map
  // front
  SfM_Data front_sfm_data;
  GoodFeatInfo gfi_front;
  if (!Load(front_sfm_data, front_sfm_data_file, ESfM_Data(VIEWS|STRUCTURE))) {
    cerr << "Failed to load " << front_sfm_data_file << endl;
    return 1;
  }
  if (!LoadData(front_sfm_data, argv[3], output_dir, "front", gfi_front)) {
    return 1;
  }
  // back
  SfM_Data back_sfm_data;
  GoodFeatInfo gfi_back;
  if (!Load(back_sfm_data, back_sfm_data_file, ESfM_Data(VIEWS|STRUCTURE))) {
    cerr << "Failed to load " << back_sfm_data_file << endl;
    return 1;
  }
  if (!LoadData(back_sfm_data, argv[3], output_dir, "back", gfi_back)) {
    return 1;
  }


  // match descs
  using MetricT = L2<unsigned char>;
  using ArrayMatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
  using Scalar = typename ArrayMatcherT::ScalarT;
  using DistanceType = typename ArrayMatcherT::DistanceType;

  ArrayMatcherT matcher;
  matcher.Build(reinterpret_cast<const Scalar*>(&gfi_back.descs[0]), gfi_back.descs.size(), 128);

  auto queries = reinterpret_cast<const Scalar*>(&gfi_front.descs[0]);
  const size_t NNN__ = 1;
  IndMatches vec_Match;
  std::vector<DistanceType> vec_Distance;
  // Search NNN__ neighbours for each query descriptor
  matcher.SearchNeighbours(queries, gfi_front.descs.size(), &vec_Match, &vec_Distance, NNN__);

  IndMatches vec_putative_match;
  vec_putative_match.reserve(vec_Match.size());
  MetricT metric;
  for (const auto &match : vec_Match) {
    if (metric(gfi_front.descs[match.i_].data(), gfi_back.descs[match.j_].data(),
      DescriptorT::static_size) < 2 * 255 * 255) {
        vec_putative_match.emplace_back(match);
    }
  }

  // for each front point, record its descs matching result
  map<int, vector<int>> point_match_info;
  if (!vec_putative_match.empty()) {
    for (const auto &putative_match : vec_putative_match) {
      int fp_index = gfi_front.points[putative_match.i_];
      int bp_index = gfi_back.points[putative_match.j_];
      point_match_info[fp_index].emplace_back(bp_index);
    }
  }

  if (point_match_info.size() == 0) {
    cerr << "Failed to find any point match" << endl;
    return 1;
  }

  using pm_info = tuple<int, int, int, int>; // fp index, bp index, score, valid descs count
  vector<pm_info> sorted_point_match_info;
  sorted_point_match_info.reserve(point_match_info.size());
  for (auto& m : point_match_info) {
    // for each front point, find its best match
    int fp_index = m.first;
    vector<int> &match_vec = m.second;
    // find the most frequent element in a vector
    sort(match_vec.begin(), match_vec.end());
    int best_score = 0;
    int best_bp_index = -1;
    for (int j = 0; j < match_vec.size(); ++j) {
      int bp_index = match_vec[j];
      int count = 1;
      for (int k = j + 1; k < match_vec.size() && match_vec[k] == bp_index; ++k) {
        ++count;
        ++j;
      }
      if (count > best_score) {
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
  for (auto& m : sorted_point_match_info) {
    auto fp_index = std::get<0>(m);
    auto bp_index = std::get<1>(m);
    auto score = std::get<2>(m);
    auto vote_num = std::get<3>(m);

    const auto &frontP = front_sfm_data.GetLandmarks().at(fp_index);
    const auto &backP = back_sfm_data.GetLandmarks().at(bp_index);
    int back_good_feat_num = 0;
    for (const auto &ob : backP.obs) {
      const auto &featmap = gfi_back.featmap_map[ob.first];
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
    for (const auto &ob : frontP.obs) {
      const auto &featmap = gfi_front.featmap_map[ob.first];
      if (featmap.find(ob.second.id_feat) != featmap.end())
        points_match_result << front_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                            << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match_result << endl;
    for (const auto &ob : backP.obs) {
      const auto &featmap = gfi_back.featmap_map[ob.first];
      if (featmap.find(ob.second.id_feat) != featmap.end())
        points_match_result << back_sfm_data.GetViews().at(ob.first).get()->s_Img_path << " "
                            << ob.second.x[0] << " " << ob.second.x[1] << " ";
    }
    points_match_result << endl;
  }
  points_match_result.close();

  return 0;
}
