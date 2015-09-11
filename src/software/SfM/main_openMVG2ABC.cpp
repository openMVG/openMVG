// Copyright (c) 2015 cpichard.
//
// // This Source Code Form is subject to the terms of the Mozilla Public
// // License, v. 2.0. If a copy of the MPL was not distributed with this
// // file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cstdlib>

// Command line
#include "third_party/cmdLine/cmdLine.h"

// OpenMVG
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/geometry/pose3.hpp"
using namespace openMVG::sfm;
using namespace openMVG::image;
using openMVG::IndexT;

// Alembic
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>
using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;
//
// The following function is taken from main_ComputeSfM_DataColor and slightly modified
//
void findTracksColors( const SfM_Data & sfm_data, 
                       const std::vector<IndexT> & trackIds,
                        std::vector<Imath::C3f> & colors)
{
    // Colorize each track
    //  Start with the most representative image
    //    and iterate to provide a color to each 3D point
    C_Progress_display my_progress_bar(sfm_data.GetLandmarks().size(),
                                       std::cout,
                                       "\nCompute scene structure color\n");

    // The track list indexes that will be colored 
    // (point removed during the process)
    std::set<IndexT> remainingTrackToColor;
    for(int i=0; i < trackIds.size(); i++)
    {
        remainingTrackToColor.insert(remainingTrackToColor.end(), i);
    }
    // Number of observation per view
    const size_t nViews = sfm_data.GetViews().size();
    std::vector<size_t> obsPerView(nViews);

    while( !remainingTrackToColor.empty() )
    {
        // Find the most representative image (for the remaining 3D points)
        //  a. Count the number of observation per view for each 3Dpoint Index
        //  b. Find the View index that is the most represented ie has the biggest
        //  number of observations
        size_t view_index = 0; // FIXME not logical, could be "no view found"
        size_t maxObs = 0;
        std::fill(obsPerView.begin(), obsPerView.end(), 0);

        for(auto iterT: remainingTrackToColor)
        {
            const size_t trackId = trackIds[iterT];
            const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
            for( Observations::const_iterator iterObs = obs.begin();
            iterObs != obs.end(); ++iterObs)
            {
                const size_t viewId = iterObs->first;
                obsPerView[viewId]++;
                if (obsPerView[viewId] >= maxObs )
                {
                    maxObs = obsPerView[viewId];
                    view_index = viewId;
                }
            }
        }

        const View * view = sfm_data.GetViews().at(view_index).get();
        const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
        Image<RGBColor> image;
        ReadImage(sView_filename.c_str(), &image);

        // Iterate through the remaining track to color
        // - look if the current view is present to color the track
        std::set<IndexT> trackProcessed;
        for(auto iterT: remainingTrackToColor)
        {
            const size_t trackId = trackIds[iterT];
            const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
            Observations::const_iterator it = obs.find(view_index);

            if(it != obs.end())
            {
                // Color the track
                const openMVG::Vec2 & pt = it->second.x;
                const RGBColor color = image(pt.y(), pt.x());
                colors[iterT] = Imath::C3f(color.r()/255.f, color.g()/255.f, color.b()/255.f);
                trackProcessed.insert(iterT);
                ++my_progress_bar;
            }
        }

        for(auto track: trackProcessed)
        {
            remainingTrackToColor.erase(track);
        }
    }
}

int main(int argc, char **argv) {
   
    // Get arguments
    CmdLine cmdLine;
    std::string sfmDataFilename;
    std::string sOutAlembic = "";

    cmdLine.add(make_option('i', sfmDataFilename, "sfmdata"));
    cmdLine.add(make_option('o', sOutAlembic, "outfile"));

    try {
        if (argc < 4) throw std::string("Invalid number of parameters in the command line.");
        cmdLine.process(argc, argv);
    } catch(const std::string &s) {
        std::cout << "openMVG to alembic\n";
        std::cout << "Usage: " << argv[0] << '\n'
        << "[-i|--sfmdata filename, the SfM_Data file to convert]\n"
        << "[-o|--outfile path]\n"
        << std::endl;
        std::cerr << s << std::endl;
        return EXIT_FAILURE;
    }

    //
    SfM_Data sfm_data;
    if (!Load(sfm_data, sfmDataFilename, ESfM_Data(ALL))) {
        std::cout << "Error: The input project file \""
                  << sfmDataFilename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }
  
    // Open Alembic archive
    // TODO: decide wether we want to keep HDF5 or switch to Ogawa 
    OArchive archive(Alembic::AbcCoreHDF5::WriteArchive(), sOutAlembic);
    OObject topObj(archive, kTop);

    //=================================================================================
    // Export points
    //=================================================================================

    std::vector<V3f> positions;
    positions.reserve(sfm_data.GetLandmarks().size());
    
    std::vector<IndexT> trackIds;
    trackIds.reserve(sfm_data.GetLandmarks().size());

    // For all the 3d points in the hash_map
    for(const auto lm : sfm_data.GetLandmarks()) {
        const openMVG::Vec3 &pt = lm.second.X;
        positions.emplace_back(pt[0], pt[1], pt[2]);
        trackIds.emplace_back(lm.first);
    }

    std::vector<Imath::C3f> colors;
    colors.resize(trackIds.size()); // R G B, 3 floats for the moment 
    findTracksColors(sfm_data, trackIds, colors);

    std::vector<Alembic::Util::uint64_t> ids(positions.size());
    std::iota(std::begin(ids), std::end(ids), 0);

    OPointsSchema::Sample psamp( std::move(V3fArraySample (positions)), std::move(UInt64ArraySample ( ids )));
    OPoints partsOut(topObj, "particleShape1");
    OPointsSchema &pSchema = partsOut.getSchema();
    pSchema.set( psamp );

    C3fArraySample val_samp( &colors[0], colors.size() );
    OC3fGeomParam::Sample color_samp(val_samp, kVertexScope);
    OCompoundProperty arbGeom = pSchema.getArbGeomParams();
    OC3fGeomParam rgbOut(arbGeom, "color", false, kVertexScope, 1);
    rgbOut.set(color_samp);
    
    //=================================================================================
    // Export cameras
    //=================================================================================
    for(const auto it: sfm_data.GetViews()) {
        const View * view = it.second.get();
        if (!sfm_data.IsPoseAndIntrinsicDefined(view))
          continue;

        // OpenMVG Camera
        const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
        auto iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
        openMVG::cameras::Pinhole_Intrinsic *cam = static_cast<openMVG::cameras::Pinhole_Intrinsic*>(iterIntrinsic->second.get());
        openMVG::Mat34 P = cam->get_projective_equivalent(pose);

        // Extract internal matrix, rotation and translation
        const openMVG::Mat3 R = pose.rotation();
        const openMVG::Vec3 center = pose.center();

        // POSE
        // Compensate translation with rotation
        // Build transform matrix
        Abc::M44d xformMatrix;
        xformMatrix[0][0]= R(0,0);
        xformMatrix[0][1]= R(0,1);
        xformMatrix[0][2]= R(0,2);
        xformMatrix[1][0]= R(1,0);
        xformMatrix[1][1]= R(1,1);
        xformMatrix[1][2]= R(1,2);
        xformMatrix[2][0]= R(2,0);
        xformMatrix[2][1]= R(2,1);
        xformMatrix[2][2]= R(2,2);
        xformMatrix[3][0]= center(0);
        xformMatrix[3][1]= center(1);
        xformMatrix[3][2]= center(2);
        xformMatrix[3][3]= 1.0;

        // Correct camera orientation for alembic
        M44d scale;
        scale[1][1] = -1;
        scale[2][2] = -1;
        
        xformMatrix = scale*xformMatrix;

        XformSample xformsample;
        xformsample.setMatrix(xformMatrix);

        stringstream ss;
        ss << stlplus::basename_part(view->s_Img_path); 
        Alembic::AbcGeom::OXform xform(topObj, "camxform_" + ss.str());
        xform.getSchema().set(xformsample);

        // Camera intrinsic parameters
        OCamera camObj(xform, "camera_" + ss.str());
        CameraSample camSample;
    
        // Take the max of the image size to handle the case where the image is in portrait mode 
        const float imgWidth = cam->w();  
        const float imgHeight = cam->h(); 
        const float sensorWidthPix = std::max(imgWidth, imgHeight);
        const float sensorHeightPix = std::min(imgWidth, imgHeight);
        const float focalLengthPix = cam->focal();
        const float dx = cam->principal_point()(0);
        const float dy = cam->principal_point()(1);
        // Use a common sensor width as we don't have this information at this point
        // We chose a full frame 24x36 camera
        const float sensorWidth = 36.0; // 36mm per default
        const float sensorHeight = sensorWidth*sensorHeightPix/sensorWidthPix;
        const float focalLength = sensorWidth*focalLengthPix/sensorWidthPix;
        // Following values are in cm, hence the 0.1 multiplier
        const float hoffset = 0.1*sensorWidth*(0.5-dx/imgWidth);
        const float voffset = 0.1*sensorHeight*(dy/imgHeight-0.5)*sensorHeightPix/sensorWidthPix;
        const float haperture = 0.1*sensorWidth*imgWidth/sensorWidthPix;
        const float vaperture = 0.1*sensorWidth*imgHeight/sensorWidthPix;

        camSample.setFocalLength(focalLength);
        camSample.setHorizontalAperture(haperture);
        camSample.setVerticalAperture(vaperture);
        camSample.setHorizontalFilmOffset(hoffset);
        camSample.setVerticalFilmOffset(voffset);

        camObj.getSchema().set( camSample );
    } 

    return EXIT_SUCCESS;
}
