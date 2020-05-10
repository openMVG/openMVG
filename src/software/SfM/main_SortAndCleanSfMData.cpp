#include <experimental/filesystem>
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"

using namespace std;
using namespace openMVG::sfm;

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    cerr << "Usage: openMVG_main_SortAndCleanSfMData sfm_data_file image_folder" << endl;
    return 1;
  }
  string sfm_data_file = argv[1];

  SfM_Data sfm_data;
  if (!Load(sfm_data, sfm_data_file, ESfM_Data(ALL)))
  {
    cerr << "Failed to load " << sfm_data_file << endl;
    return 1;
  }

  std::vector<std::string> view_names;
  for (auto &p : std::experimental::filesystem::directory_iterator(argv[2]))
  {
    if (std::experimental::filesystem::is_regular_file(p))
    {
      view_names.push_back(p.path().filename());
    }
  }

  SortAndCleanSfMData(sfm_data, view_names);

  Save(sfm_data, sfm_data_file, ESfM_Data(ALL));

  return 0;
}
