#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

int main(int argc, const char *argv[]) {
    if (argc < 3) {
        std::cout << "Usage: undistort sfm_data_json input_output_list" << std::endl;
        return 1;
    }

    SfM_Data sfm_data;
    if (!Load(sfm_data, argv[1], ESfM_Data(INTRINSICS))) {
        std::cerr << std::endl
            << "The input SfM_Data file \"" << argv[1] << "\" cannot be read." << std::endl;
        return 1;
    }

    if (sfm_data.intrinsics.size() != 1) {
        std::cerr << std::endl
            << "SfM_Data can only contain 1 intrinsic" << std::endl;
        return 1;
    }

    auto intrinsic = sfm_data.intrinsics.begin()->second;
    cv::Mat map1(intrinsic->h_, intrinsic->w_, CV_32F);
    cv::Mat map2 = map1.clone();
    for (int i = 0; i < map1.rows; ++i) {
        for (int j = 0; j < map2.cols; ++j) {
            const Vec2 undisto_pix(j, i);
            const Vec2 disto_pix = intrinsic->get_d_pixel(undisto_pix);
            map1.at<float>(i, j) = (float)disto_pix[0];
            map2.at<float>(i, j) = (float)disto_pix[1];
        }
    }

    std::ifstream list(argv[2]);
    if (!list.good()) {
        std::cerr << std::endl
            << "The input list file \"" << argv[2] << "\" cannot be read." << std::endl;
        return 1;
    }

    std::string input, output;
    while (list >> input && list >> output) {
        cv::Mat img = cv::imread(input);
        cv::Mat undistorted;
        cv::remap(img, undistorted, map1, map2, cv::INTER_LANCZOS4);
        cv::imwrite(output, undistorted);
    }

    return 0;
}
