#ifndef CSIFT_IMAGE_DESCRIBER_HPP
#define CSIFT_IMAGE_DESCRIBER_HPP

#include <tuple>
#include <vector>
#include "openMVG/image/image.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/feature.hpp"
#include "third_party/cudasift/cudaImage.h"
#include "third_party/cudasift/cudaSift.h"
#include "third_party/cudasift/csift_keypoint.hpp"
#include <iostream>
#include <numeric>
#include <cereal/cereal.hpp>
#include <Eigen/Core>

namespace openMVG {
namespace features {

class CSIFT_Image_describer : public Image_describer
{
public:
    struct Params
    {
        Params(
        float initBlur = 1.0f,
        float thresh = 3.5f
        ):
        initBlur_(initBlur),
        thresh_(thresh) {}

        template<class Archive>
        void serialize(Archive & ar)
        {
            ar(
                cereal::make_nvp("initBlur", initBlur_),
                cereal::make_nvp("thresh", thresh_));

        }
        float initBlur_;
        float thresh_;
    };

    CSIFT_Image_describer
    (
        const Params params = Params()
    )
    :Image_describer(), params_(params)
    {}

    bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
    {
        switch(preset)
        {
          case NORMAL_PRESET:
            params_.thresh_ = 2.5f;
          break;
          case HIGH_PRESET:
            params_.thresh_ = 3.5f;
          break;
          case ULTRA_PRESET:
            params_.thresh_ = 4.00f;
          break;
          default:
            return false;
        }
        return true;
    }


    bool Describe(
        const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr) override
    {
        InitCuda(0);
        CudaImage img1;
        const int w = image.Width(), h = image.Height();


        // Convert to float in range [0;1]
        const image::Image<float> If(image.GetMat().cast<float>());
        img1.Allocate(w, h, iAlignUp(w, 128), false, NULL, (float*)If.data());
        img1.Download();
        Allocate(regions);
        CSIFT_Regions * regionsCasted = dynamic_cast<CSIFT_Regions*>(regions.get());
        {
            SiftData siftData1;
            InitSiftData(siftData1, 32768, true, true);
	          //std::string filename = "/home/sai/lambe.png";
	    //image::WriteImage(filename.c_str(), image.data());
	    std::cout << "Config is" << params_.initBlur_ << "and " << params_.thresh_ << std::endl;
            ExtractSift(siftData1, img1, 5, params_.initBlur_, params_.thresh_, 0.0f, false);
            int numPts1 = siftData1.numPts;
	    std::cout << "Number of features detected: " << numPts1 << std::endl;
	    //std::cout << "Image data" << std::endl;
	    //std::cout << image.data();
            SiftPoint *sift1 = siftData1.h_data;
            std::vector< csift::Keypoint > keys;
            Descriptor<unsigned char, 128> descriptor;
            for (int i=0;i<numPts1;i++)
            {
                csift::Keypoint k;
                k.x = sift1[i].xpos;
                k.y = sift1[i].ypos;
                k.sigma = sift1[i].scale;
                k.theta = sift1[i].orientation;
                k.descr = Eigen::Map<Eigen::Matrix<float,128,1>> (sift1[i].data);
                keys.push_back(k);
            }
            std::vector< csift::Keypoint > keypoints;
            std::move(keys.begin(), keys.end(), std::back_inserter(keypoints));

            for (const auto & k : keypoints)
            {
                Descriptor<unsigned char, 128> descriptor;
                descriptor << (k.descr.cast<unsigned char>());
                {
                  for (int ctr=0;ctr<128;++ctr)
                    descriptor[ctr] = static_cast<unsigned char>(512.f*k.descr[ctr]);


                    regionsCasted->Descriptors().emplace_back(descriptor);
                    regionsCasted->Features().emplace_back(k.x, k.y, k.sigma, k.theta);
                }
            }
        }
        return true;
    }

    void Allocate(std::unique_ptr<Regions> &regions) const override
    {
        regions.reset( new CSIFT_Regions );
    }

    template<class Archive>
    void serialize( Archive & ar )
    {
        ar(cereal::make_nvp("params", params_));
    }


private:
    Params params_;
};

}
}

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::CSIFT_Image_describer, "CSIFT_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::CSIFT_Image_describer)
#endif // OPENMVG_FEATURES_SIFT_SIFT_ANATOMY_IMAGE_DESCRIBER_HPP
