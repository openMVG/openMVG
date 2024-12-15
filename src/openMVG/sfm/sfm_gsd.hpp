
#include <vector>

#include "openMVG/sfm/sfm_data.hpp"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

namespace openMVG {
namespace sfm {

bool AreSame(double a, double b);

Vec2 get_ground_sampling_distance(
    const SfM_Data & sfm_data, 
    const Vec2 & sensor_size
);

std::vector<Vec2> getRealGroundCoordinates(
    IndexT view_id, 
    const SfM_Data & data,
    const double & average_height
);

}
}