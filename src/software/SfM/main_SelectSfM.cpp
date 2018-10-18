#include <fstream>
#include <unordered_map>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

using namespace std;
using namespace openMVG;
using namespace openMVG::sfm;

struct InputSlice {
  enum SliceType { kAll, kSlice, kList };

  string sfm_json;
  SliceType type{kAll};
  int start{0}, end{0};
  vector<string> filenames;

  bool Contains(const shared_ptr<View> &view, size_t size) {
    switch (type) {
    case kAll:
      return true;
    case kSlice: {
      int true_start = start < 0 ? start + (int)size : start;
      int true_end = end <= 0 ? end + (int)size : end;
      int index = (int)view->id_view;
      return index >= true_start && index < true_end;
    }
    case kList:
      return find(filenames.begin(), filenames.end(), view->s_Img_path) != filenames.end();
    default:
      break;
    }
    return false;
  }
};

bool ParseInput(int argc, char *argv[], int &i, InputSlice &input) {
  input.sfm_json = argv[i];
  input.type = InputSlice::kAll;
  ++i;

  while (i < argc) {
    string argi = argv[i];
    if (argi[0] != '-')
      return true;

    if (argi == "-s") {
      ++i;
      if (i >= argc)
        return false;
      input.type = InputSlice::kSlice;
      input.start = stoi(argv[i]);
      ++i;
    } else if (argi == "-e") {
      ++i;
      if (i >= argc)
        return false;
      input.type = InputSlice::kSlice;
      input.end = stoi(argv[i]);
      ++i;
    } else if (argi == "-l") {
      ++i;
      if (i >= argc)
        return false;
      ifstream f(argv[i]);
      if (!f.good())
        return false;
      input.type = InputSlice::kList;
      string line;
      while (getline(f, line))
        input.filenames.push_back(std::move(line));
      ++i;
    } else {
      return false;
    }
  }
  return true;
}

bool ParseCmdLine(int argc, char *argv[], string &output_filename, string &root_path, vector<InputSlice> &inputs) {
  if (argc < 3)
    return false;
  output_filename = argv[1];
  root_path = argv[2];

  inputs.clear();
  int i = 3;
  while (i < argc) {
    InputSlice input;
    if (!ParseInput(argc, argv, i, input))
      return false;
    inputs.push_back(input);
  }
  return true;
}

int main(int argc, char *argv[]) {
  string output_filename, root_path;
  vector<InputSlice> inputs;
  if (!ParseCmdLine(argc, argv, output_filename, root_path, inputs)) {
    cerr << "Usage: openMVG_main_SelectSfM output_filename root_path "
      "input_filename1 [-s start_index -e end_index] [-l list_filename] input_filename2 ..." << endl;
    return 1;
  }

  vector<tuple<string, int, int>> stats;

  SfM_Data output_sfm;
  output_sfm.s_root_path = root_path;
  IndexT view_index = 0, intrinsics_index = 0;
  for (auto &input : inputs) {
    SfM_Data input_sfm;
    if (!Load(input_sfm, input.sfm_json, ESfM_Data(VIEWS|INTRINSICS))) {
      cerr << "Failed to load " << input.sfm_json << endl;
      continue;
    }

    unordered_map<IndexT, IndexT> intrinsics_id_map;
    int view_count = 0, intrinsics_count = 0;

    for (auto view_it = input_sfm.views.begin(); view_it != input_sfm.views.end(); ++view_it) {
      if (!input.Contains(view_it->second, input_sfm.views.size()))
        continue;

      IndexT ori_intrin_index = view_it->second->id_intrinsic;
      auto intrin_it = input_sfm.intrinsics.find(ori_intrin_index);
      if (intrin_it != input_sfm.intrinsics.end()) {
        auto map_it = intrinsics_id_map.find(ori_intrin_index);
        if (map_it == intrinsics_id_map.end()) {
          // copy intrinsics to output sfm, and record their new index
          output_sfm.intrinsics[intrinsics_index] = intrin_it->second;
          intrinsics_id_map[ori_intrin_index] = intrinsics_index;
          ++intrinsics_index;
          ++intrinsics_count;
        }

        // copy view to output sfm, change its ids
        output_sfm.views[view_index] = view_it->second;
        view_it->second->id_view = view_index;
        view_it->second->id_pose = view_index;
        view_it->second->id_intrinsic = intrinsics_id_map[ori_intrin_index];
        ++view_index;
        ++view_count;
      } else {
        cerr << "Discarding view " << view_it->second->s_Img_path << " in " << input.sfm_json << endl;
      }
    }

    stats.emplace_back(input.sfm_json, view_count, intrinsics_count);
  }

  if (!Save(output_sfm, output_filename, ESfM_Data(VIEWS|INTRINSICS))) {
    cerr << "Failed to save to " << output_filename << std::endl;
    return 1;
  }

  cout << endl << "SelectSfM report:" << endl;
  for (auto &stat : stats) {
    cout << "selected #" << get<1>(stat) << " views and #" << get<2>(stat) << " intrinsics from " << get<0>(stat) << endl;
  }
  cout << "total #" << view_index << " views and #" << intrinsics_index << " intrinsics" << endl;

  return 0;
}
