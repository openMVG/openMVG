#include <cereal/archives/json.hpp>

#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/features/image_describer_akaze_io.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer_io.hpp"
#include "nonFree/sift/SIFT_describer_io.hpp"
#include "nonFree/rich_sift/Rich_SIFT_describer_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_io.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace std;
using namespace stlplus;
using namespace openMVG::sfm;
using namespace openMVG::features;
using namespace openMVG::image;

bool find_feature_files(const string &image_name, const vector<string> &dirs,
  string &feat_file, string &desc_file) {
  string basename = basename_part(image_name);
  for (auto &dir : dirs) {
    feat_file = create_filespec(dir, basename, "feat");
    desc_file = create_filespec(dir, basename, "desc");
    if (file_exists(feat_file) && file_exists(desc_file))
      return true;
  }
  return false;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: openMVG_main_MaskFeatures input_sfm_filename output_dir feature_dir1 feature_dir2 ..." << endl;
    return 1;
  }

  const char *input_sfm = argv[1];
  SfM_Data sfm_data;
  if (!Load(sfm_data, input_sfm, VIEWS)) {
    cerr << "Failed to load " << input_sfm << endl;
    return 1;
  }
  string output_dir = argv[2];

  vector<string> feature_dirs;
  for (int i = 3; i < argc; ++i)
    feature_dirs.push_back(argv[i]);

  unique_ptr<Image_describer> describer;
  unique_ptr<Regions> regions;
  for (auto it = sfm_data.views.begin(); it != sfm_data.views.end(); ++it) {
    string feat_file, desc_file;
    if (!find_feature_files(it->second->s_Img_path, feature_dirs, feat_file, desc_file)) {
      cerr << "Failed to find feature files for " << it->second->s_Img_path << endl;
      continue;
    }

    // check describer type
    string folder = folder_part(feat_file);
    string describer_file = create_filespec(folder, "image_describer", "json");
    ifstream istream(describer_file);
    if (!istream.good()) {
      cerr << "Missing image describer " << describer_file << endl;
      continue;
    }

    describer.reset();
    try {
      cereal::JSONInputArchive archive(istream);
      archive(cereal::make_nvp("image_describer", describer));
    } catch (const cereal::Exception &e) {
      cerr << e.what() << endl << "Failed to load image describer" << endl;
      continue;
    }

    if (regions) {
      auto this_regions = describer->Allocate();
      if (this_regions->Type_id() != regions->Type_id() ||
        this_regions->DescriptorLength() != regions->DescriptorLength()) {
        cerr << "Different describer type, abandoning " << folder << endl;
        continue;
      }
    } else {
      regions = describer->Allocate();
    }

    // load features
    if (!regions->Load(feat_file, desc_file)) {
      cerr << "Failed to load feature " << feat_file << endl;
      continue;
    }

    // load mask
    bool has_mask = true;
    string mask_file = create_filespec(sfm_data.s_root_path, it->second->s_Img_path, ".mask.png");
    Image<unsigned char> mask_image;
    if (!ReadImage(mask_file.c_str(), &mask_image)) {
      cerr << "\nFailed to read mask " << mask_file << endl;
      has_mask = false;
    }

    // copy features in the mask to new regions
    auto new_regions = describer->Allocate();
    vector<size_t> old_indices;
    for (size_t i = 0; i < regions->RegionCount(); ++i) {
      auto pos = regions->GetRegionPosition(i);
      if (has_mask && mask_image((int)pos.y(), (int)pos.x())) {
        regions->CopyRegion(i, new_regions.get());
        old_indices.push_back(i);
      }
    }

    describer->Describe(mask_image, new_regions);


    // save features
    string basename = basename_part(it->second->s_Img_path);
    feat_file = create_filespec(output_dir, basename, "feat");
    desc_file = create_filespec(output_dir, basename, "desc");
    if (!describer->Save(new_regions.get(), feat_file, desc_file)) {
      cerr << "Failed to save features to " << feat_file << endl;
      continue;
    }

    // save index map
    string map_file = create_filespec(output_dir, basename, "featmap");
    ofstream ostream(map_file);
    if (!ostream.good()) {
      cerr << "Failed to open map file " << map_file << endl;
      continue;
    }
    for (auto i : old_indices)
      ostream << i << endl;

    // print stats
    cout << "Selected #" << new_regions->RegionCount() << " in #"
      << regions->RegionCount() << " features from " << it->second->s_Img_path << endl;
  }

  if (describer) {
    string describer_file = create_filespec(output_dir, "image_describer", "json");
    ofstream stream(describer_file);
    if (!stream.good()) {
      cerr << "Failed to open describer file " << describer_file << endl;
      return 1;
    }

    {
      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", describer));
      regions = describer->Allocate();
      archive(cereal::make_nvp("regions_type", regions));
    }

    cout << "Describer saved to " << describer_file << endl;
    return 0;
  }

  cerr << "No valid describer" << endl;
  return 1;
}
