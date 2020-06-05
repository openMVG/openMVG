#include "openMVG/spherical/spherical.hpp"
#include "openMVG/image/image_io.hpp"

#include "testing/testing.h"

#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::spherical;

TEST(Spherical, EquirectToTangent) {
  Image<RGBColor> image;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  std::cout << png_filename << std::endl;
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  std::vector<Image<RGBColor>> tangent_images;
  CreateTangentImages(image, 9L, 0L, tangent_images);

  for (size_t i = 0; i < tangent_images.size(); i++) {
    const std::string out_filename =
        ("test_tangent images_" + std::to_string(i) + ".png");
    EXPECT_TRUE(WriteImage(out_filename.c_str(), tangent_images[i]));
  }
}

TEST(Spherical, TangentToEquirect) {
  Image<RGBColor> image;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  std::cout << png_filename << std::endl;
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  std::vector<Image<RGBColor>> tangent_images;
  CreateTangentImages(image, 9L, 0L, tangent_images);

  // TangentToEquirectangular(tangent_images, )
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */