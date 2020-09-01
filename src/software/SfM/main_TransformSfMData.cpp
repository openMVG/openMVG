#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"

using namespace std;
using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char *argv[])
{
  if (argc != 19)
  {
    cerr << "Usage: " << argv[0] << " input_sfm_data_file output_sfm_data_file transform_matrix" << endl;
    return 1;
  }

  std::string input_sfm_data_file = argv[1];
  std::string output_sfm_data_file = argv[2];

  SfM_Data sfm_data;
  if (!Load(sfm_data, input_sfm_data_file, ESfM_Data(ALL)))
  {
    cerr << "Faied to load " << input_sfm_data_file << endl;
    return 1;
  }

  Mat4 transform;
  for (int row = 0; row < 4; ++row)
  {
    for (int col = 0; col < 4; ++col)
    {
      transform(row, col) = stod(argv[row * 4 + col + 3]);
    }
  }
  TransformSfMData(transform, sfm_data);

  if (!Save(sfm_data, output_sfm_data_file, ESfM_Data(ALL)))
  {
    cerr << "Faied to save sfm data to " << output_sfm_data_file << endl;
    return 1;
  }

  return 0;
}
