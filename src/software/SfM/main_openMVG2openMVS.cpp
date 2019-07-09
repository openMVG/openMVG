// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 cDc <cdc.seacave@gmail.com>, Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#define _USE_EIGEN
#include "InterfaceMVS.h"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress_display.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include <cstdlib>
#include <ctype.h>
#include <string>
#include <array>
#include <unordered_map>
#include <memory>

template<typename Image, typename RemapType, typename SamplerType>
void UndistortImage(
    const Image& imageIn,
    Image& imageUd,
    const RemapType& mapRow,
    const RemapType& mapCol,
    const SamplerType& sampler,
    typename Image::Tpixel fillcolor = typename Image::Tpixel(0))
{
    int width = imageIn.Width();
    int height = imageIn.Height();

    imageUd.resize(width, height, true, fillcolor);

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for (int j = 0; j < height; ++j)
        for (int i = 0; i < width; ++i)
        {
            // compute coordinates with distortion
            double row = (double)mapRow(j, i);
            double col = (double)mapCol(j, i);
            // pick pixel if it is in the image domain
            if (imageIn.Contains(row, col))
            {
                imageUd(j, i) = sampler(imageIn, row, col);
            }
        }
}

template<typename RemapType>
void CalcCameraDistortionMap(
    const openMVG::cameras::IntrinsicBase* cam,
    RemapType& outMapRow,
    RemapType& outMapCol
)
{
    int width = cam->w();
    int height = cam->h();

    outMapRow.resize(width, height);
    outMapCol.resize(width, height);

#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const Vec2 undistoPix(i, j);
            // compute coordinates with distortion
            const Vec2 distoPix = cam->get_d_pixel(undistoPix);

            outMapRow(j, i) = distoPix[1];
            outMapCol(j, i) = distoPix[0];
        }
}

bool exportToUndistortDepths(
    const std::string& sInDir,
    const SfM_Data& sfmData,
    const std::string& sOutDir
)
{
    if (sfmData.GetViews().empty())
        return false;
    int viewIdx = -1;
    for (const auto &v : sfmData.GetViews())
    {
        if (sfmData.IsPoseAndIntrinsicDefined(v.second.get()))
        {
            viewIdx = v.first;
            break;
        }
    }
    if (viewIdx == -1)
        return false;
    const auto& view = sfmData.GetViews().at(viewIdx);
    const auto& ext = stlplus::extension_part(view->s_Img_path);

    const auto* cam = sfmData.GetIntrinsics().at(view->id_intrinsic).get();
    if (cam == nullptr)
        return false;

    if (!stlplus::is_folder(sOutDir))
    {
        stlplus::folder_create(sOutDir);
        if (!stlplus::is_folder(sOutDir))
        {
            std::cerr << "Cannot access to one of the desired output directory" << std::endl;
            return false;
        }
    }

    std::unordered_map<const openMVG::cameras::IntrinsicBase*, std::array<Image<double>, 2>> mapCamRemap;

    bool ret = true;

    std::vector<std::string> vecFile = stlplus::folder_files(sInDir);
    std::sort(vecFile.begin(), vecFile.end());
    for (const auto& file : vecFile)
    {
        if (stlplus::extension_part(file) != "pfm")
            continue;

        const std::string srcFile = stlplus::create_filespec(sInDir, file);
        if (!stlplus::is_file(srcFile))
        {
            std::cout << "Cannot read the corresponding depth file: " << srcFile << std::endl;
            ret = false;
            break;
        }

        std::string outputFile = stlplus::create_filespec(sOutDir, file);

        // export undistorted depths
        if (cam->have_disto())
        {
            auto iterMapCam = mapCamRemap.find(cam);
            if (mapCamRemap.end() == iterMapCam)
            {
                // map row and col have not been cached
                CalcCameraDistortionMap(cam, mapCamRemap[cam][0], mapCamRemap[cam][1]);
            }

            auto& mapRow = mapCamRemap[cam][0];
            auto& mapCol = mapCamRemap[cam][1];

            // undistort depth and save it
            Image<float> depthFloat, depthFloatUd;
            ReadImage(srcFile.c_str(), &depthFloat);
            UndistortImage(depthFloat, depthFloatUd, mapRow, mapCol, image::Sampler2d<image::SamplerNearest>(), 0.0f);
            WriteImage(outputFile.c_str(), depthFloatUd);
        }
        else
        {
            // just copy image
            stlplus::file_copy(srcFile, outputFile);
        }
    }

    return ret;
}

bool exportToUndistortImages(
    const SfM_Data& sfmData,
    const std::string& sOutDir
)
{
    if (!stlplus::is_folder(sOutDir))
    {
        stlplus::folder_create(sOutDir);
        if (!stlplus::is_folder(sOutDir))
        {
            std::cerr << "Cannot access to one of the desired output directory" << std::endl;
            return false;
        }
    }

    std::unordered_map<const openMVG::cameras::IntrinsicBase*, std::array<Image<double>, 2>> mapCamRemap;

    bool ret = true;

    for (const auto& view : sfmData.GetViews())
    {
        const std::string srcImage = stlplus::create_filespec(sfmData.s_root_path, view.second->s_Img_path);
        const std::string ext = stlplus::extension_part(view.second->s_Img_path);

        const std::string outputImage = stlplus::create_filespec(sOutDir, view.second->s_Img_path);

        if (!stlplus::is_file(srcImage))
        {
            std::cout << "Cannot read the corresponding image: " << srcImage << std::endl;
            ret = false;
            break;
        }
        if (sfmData.IsPoseAndIntrinsicDefined(view.second.get()))
        {
            // export undistorted images
            const openMVG::cameras::IntrinsicBase* cam = sfmData.GetIntrinsics().at(view.second->id_intrinsic).get();

            if (cam->have_disto())
            {
                auto iterMapCam = mapCamRemap.find(cam);
                if (mapCamRemap.end() == iterMapCam)
                {
                    // map row and col have not been cached
                    CalcCameraDistortionMap(cam, mapCamRemap[cam][0], mapCamRemap[cam][1]);
                }

                auto& mapRow = mapCamRemap[cam][0];
                auto& mapCol = mapCamRemap[cam][1];

                // undistort image and save it
                Image<openMVG::image::RGBColor> imageRGB, imageRGBUd;
                ReadImage(srcImage.c_str(), &imageRGB);
                UndistortImage(imageRGB, imageRGBUd, mapRow, mapCol, image::Sampler2d<image::SamplerLinear>(), BLACK);
                WriteImage(outputImage.c_str(), imageRGBUd);
                // undistort mask and save it
                auto maskFileName = srcImage + ".mask.png";
                if (stlplus::is_file(maskFileName)) {
                    Image<openMVG::image::RGBColor> maskRGB, maskRGBUd;
                    ReadImage(maskFileName.c_str(), &maskRGB);
                    UndistortImage(maskRGB, maskRGBUd, mapRow, mapCol, image::Sampler2d<image::SamplerNearest>(), BLACK);
                    WriteImage((outputImage + ".mask.png").c_str(), maskRGBUd);
                }
                // undistort depth and save it
                auto depthFileBasename = view.second->s_Img_path;
                auto pos = depthFileBasename.find("rgb");
                if (pos != std::string::npos) {
                    depthFileBasename.replace(pos, 3, "depth");
                }
                pos = depthFileBasename.rfind(std::string(".") + ext);
                if (pos != std::string::npos) {
                    depthFileBasename.replace(pos, 4, ".pfm");
                }
                auto depthFileName = stlplus::create_filespec(sfmData.s_root_path, depthFileBasename);

                if (depthFileBasename != view.second->s_Img_path && stlplus::is_file(depthFileName)) {
                    Image<float> depthFloat, depthFloatUd;
                    ReadImage(depthFileName.c_str(), &depthFloat);
                    UndistortImage(depthFloat, depthFloatUd, mapRow, mapCol, image::Sampler2d<image::SamplerNearest>(), 0.0f);
                    WriteImage((outputImage + ".ref.pfm").c_str(), depthFloatUd);
                }
            }
            else
            {
                // just copy image
                stlplus::file_copy(srcImage, outputImage);
                stlplus::file_copy(srcImage + ".mask.png", outputImage + ".mask.png");
                stlplus::file_copy(srcImage + ".ref.pfm", outputImage + ".ref.pfm");
            }
        }
        else
        {
            // just copy the image
            stlplus::file_copy(srcImage, outputImage);
        }
    }

    return ret;
}

bool exportToOpenMVS(
    const SfM_Data& sfmData,
    const std::string& sOutFile
)
{
    // Export data :
    MVS::Interface scene;
    size_t nPoses(0);
    const uint32_t nViews((uint32_t)sfmData.GetViews().size());

    C_Progress_display myProgressBar(nViews);

    // OpenMVG can have not contiguous index, use a map to create the required OpenMVS contiguous ID index
    std::map<openMVG::IndexT, uint32_t> mapIntrinsic, mapView;

    // define a platform with all the intrinsic group
    for (const auto& intrinsic : sfmData.GetIntrinsics())
    {
        if (isPinhole(intrinsic.second->getType()))
        {
            const Pinhole_Intrinsic* cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());
            if (mapIntrinsic.count(intrinsic.first) == 0)
                mapIntrinsic.insert(std::make_pair(intrinsic.first, scene.platforms.size()));
            MVS::Interface::Platform platform;
            // add the camera
            MVS::Interface::Platform::Camera camera;
            camera.K = cam->K();
            // sub-pose
            camera.R = Mat3::Identity();
            camera.C = Vec3::Zero();
            platform.cameras.push_back(camera);
            scene.platforms.push_back(platform);
        }
    }

    // define images & poses
    scene.images.reserve(nViews);
    for (const auto& view : sfmData.GetViews())
    {
        mapView[view.first] = scene.images.size();
        MVS::Interface::Image image;
        // store relative path because we going to feed mvs resized images
        image.name = view.second->s_Img_path;
        image.platformID = mapIntrinsic.at(view.second->id_intrinsic);
        MVS::Interface::Platform& platform = scene.platforms[image.platformID];
        image.cameraID = 0;

        if (sfmData.IsPoseAndIntrinsicDefined(view.second.get()))
        {
            MVS::Interface::Platform::Pose pose;
            image.poseID = platform.poses.size();
            const openMVG::geometry::Pose3 poseMVG(sfmData.GetPoseOrDie(view.second.get()));
            pose.R = poseMVG.rotation();
            pose.C = poseMVG.center();
            platform.poses.push_back(pose);
            ++nPoses;
        }
        else
        {
            // image have not valid pose, so set an undefined pose
            image.poseID = NO_ID;
        }
        scene.images.emplace_back(image);
        ++myProgressBar;
    }

    // define structure
    scene.vertices.reserve(sfmData.GetLandmarks().size());
    for (const auto& vertex : sfmData.GetLandmarks())
    {
        const Landmark& landmark = vertex.second;
        MVS::Interface::Vertex vert;
        MVS::Interface::Vertex::ViewArr& views = vert.views;
        for (const auto& observation : landmark.obs)
        {
            const auto it(mapView.find(observation.first));
            if (it != mapView.end()) {
                MVS::Interface::Vertex::View view;
                view.imageID = it->second;
                view.confidence = 0;
                views.push_back(view);
            }
        }
        if (views.size() < 2)
            continue;
        std::sort(
            views.begin(), views.end(),
            [](const MVS::Interface::Vertex::View & view0, const MVS::Interface::Vertex::View & view1)
            {
                return view0.imageID < view1.imageID;
            }
        );
        vert.X = landmark.X.cast<float>();
        scene.vertices.push_back(vert);
    }

    // normalize camera intrinsics
    for (size_t p = 0; p < scene.platforms.size(); ++p)
    {
        MVS::Interface::Platform& platform = scene.platforms[p];
        for (size_t c = 0; c < platform.cameras.size(); ++c) {
            MVS::Interface::Platform::Camera& camera = platform.cameras[c];
            // find one image using this camera
            MVS::Interface::Image* pImage(nullptr);
            for (MVS::Interface::Image& image : scene.images)
            {
                if (image.platformID == p && image.cameraID == c && image.poseID != NO_ID)
                {
                    pImage = &image;
                    break;
                }
            }
            if (pImage == nullptr)
            {
                std::cerr << "error: no image using camera " << c << " of platform " << p << std::endl;
                continue;
            }
            // read image meta-data
            ImageHeader imageHeader;
            ReadImageHeader(stlplus::create_filespec(sfmData.s_root_path, pImage->name).c_str(), &imageHeader);
            const double fScale(1.0 / std::max(imageHeader.width, imageHeader.height));
            camera.K(0, 0) *= fScale;
            camera.K(1, 1) *= fScale;
            camera.K(0, 2) *= fScale;
            camera.K(1, 2) *= fScale;
        }
    }

    // write OpenMVS data
    if (!MVS::ARCHIVE::SerializeSave(scene, sOutFile))
        return false;

    std::cout
        << "Scene saved to OpenMVS interface format:\n"
        << " #platforms: " << scene.platforms.size() << std::endl;
    for (int i = 0; i < scene.platforms.size(); ++i)
    {
        std::cout << "  platform ( " << i << " ) #cameras: " << scene.platforms[i].cameras.size() << std::endl;
    }
    std::cout
        << "  " << scene.images.size() << " images (" << nPoses << " calibrated)\n"
        << "  " << scene.vertices.size() << " Landmarks\n";
    return true;
}


int main(int argc, char *argv[])
{
    CmdLine cmd;
    std::string sSfMDataFilename;
    std::string sInDir = "";
    std::string sOutFile = "";
    std::string sOutDir = "";

    cmd.add( make_option('i', sSfMDataFilename, "sfmdata") );
    cmd.add( make_option('I', sInDir, "indir"));
    cmd.add( make_option('o', sOutFile, "outfile") );
    cmd.add( make_option('d', sOutDir, "outdir") );

    try {
        if (argc == 1) throw std::string("Invalid command line parameter.");
        cmd.process(argc, argv);
    } catch (const std::string& s) {
        std::cerr << "Usage: " << argv[0] << '\n'
        << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
        << "[-I]--indir] depth images path to convert\n"
        << "[-o|--outfile] OpenMVS scene file\n"
        << "[-d|--outdir] undistorted images path\n"
        << std::endl;

        std::cerr << s << std::endl;
        return EXIT_FAILURE;
    }

    if ((0 == sOutFile.size()) && (0 == sOutDir.size()))
    {
        std::cerr << std::endl
            << "At least one of OpenMVS scene file and undistorted images path should be assigned"
            << std::endl;
    }

    if ((sOutFile.size() > 0) && (stlplus::extension_part(sOutFile) != "mvs")) {
        std::cerr << std::endl
        << "Invalid output file extension: " << sOutFile << std::endl
        << "You must use a filename with a .mvs extension." << std::endl;
        return EXIT_FAILURE;
    }

    // Read the input SfM scene
    SfM_Data sfmData;
    if (!Load(sfmData, sSfMDataFilename, ESfM_Data(ALL))) {
        std::cerr << std::endl
        << "The input SfM_Data file \""<< sSfMDataFilename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }

    bool onlyExportDepth = !sInDir.empty();
    bool exportOpenMVS = sInDir.empty();

    if (sOutDir.size() > 0)
    {
        if (onlyExportDepth)
        {
            if (!exportToUndistortDepths(sInDir, sfmData, sOutDir))
            {
                std::cerr << std::endl
                    << "Error occurred during generation of undistorted images." << std::endl;
                return EXIT_FAILURE;
            }
        }
        else
        {
            if (!exportToUndistortImages(sfmData, sOutDir))
            {
                std::cerr << std::endl
                    << "Error occurred during generation of undistorted images." << std::endl;
                return EXIT_FAILURE;
            }
        }
    }

    if (exportOpenMVS && sOutFile.size() > 0)
    {
        if (!exportToOpenMVS(sfmData, sOutFile))
        {
            std::cerr << std::endl
                << "The output openMVS scene file can not be written." << std::endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
