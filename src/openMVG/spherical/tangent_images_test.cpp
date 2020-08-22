#include "openMVG/spherical/tangent_images.hpp"
#include "nonFree/sift/SIFT_describer.hpp"
#include "openMVG/image/image_io.hpp"

#include "testing/testing.h"

#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::spherical;
using namespace openMVG::features;

TEST(Spherical, EquirectToTangent) {
  // Load the test equirectangular image
  Image<RGBColor> image;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  // Instantiate a TangentImage object that defines the relationship between the dimension of the equirectangular image and those of the tangent images we want to create
  TangentImages tangent_images(0, 9, image.Height(), image.Width());

  // For fun, this is the FOV of each tangent image we're going to create
  std::cout << "FOV: " << tangent_images.FOV() << " degrees" << std::endl;

  // Create the tangent images
  std::vector<Image<RGBColor>> t_images;
  tangent_images.CreateTangentImages(image, t_images);

  // Now write the create tangent images to file
  for (size_t i = 0; i < t_images.size(); i++) {
    const std::string out_filename =
        ("../openMVG/test_tangent images_" + std::to_string(i) + ".png");
    EXPECT_TRUE(WriteImage(out_filename.c_str(), t_images[i]));
  }
}


TEST(Spherical, DescribeTangentImages) {
  // Load the test equirectangular image as a grayscale image
  Image<unsigned char> image;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  // Instantiate the TangentImages object
  TangentImages tangent_images(0, 9, image.Height(), image.Width());

  // Create the tangent images
  std::vector<Image<unsigned char>> t_images;
  tangent_images.CreateTangentImages(image, t_images);

  // Instantiate the image describer of choice (SIFT used here)
  std::unique_ptr<Image_describer> image_describer;
  image_describer.reset(
      new SIFT_Image_describer(SIFT_Image_describer::Params(), true));

  // Also create a storage vector for the regions of each tangent image
  std::vector<std::unique_ptr<Regions>> all_regions;

  // Go through each tangent image and run the feature detector
  for (size_t i = 0; i < t_images.size(); i++){
    std::unique_ptr<Regions> regions = image_describer->Describe(t_images[i]);

    // We can save the individual outputs for each tangent image if we need it
    const std::string sFeat = "../openMVG/image" + std::to_string(i) + ".feat";
    const std::string sDesc = "../openMVG/image" + std::to_string(i) + ".desc";
    image_describer->Save(regions.get(), sFeat, sDesc);

    // Store the features for the set of tangent images in a vector. Note that, because we use std::unique_ptr to reference the features, we must use C++11's move semantics to store it in the vector
    all_regions.push_back(std::move(regions));
  }

  /* At this point, we have points and descriptors computed, but the coordinates all correspond to UV coordinates on the tangent image. Most likely, we want to convert those coordinates back to the equirectangular image. Additionally, we have Regions objects for each tangent image, but we really just want one Regions object that contains all the features for the equirectangular image. This next step will do these desired conversions. */
  
  // std::unique_ptr<Regions> equirect_regions =
      // tangent_images.ConvertTangentImageFeaturesToEquirectangular(all_regions);
}



/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */