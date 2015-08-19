#include <cstdlib>

// Command line
#include "third_party/cmdLine/cmdLine.h"

// OpenMVG
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/geometry/pose3.hpp"
using namespace openMVG::sfm;

// Alembic
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>
using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;


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
        std::cout << "openMVG to alembic" << std::endl;
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
        cout << "Error: The input project file \""<< sfmDataFilename << "\" cannot be read." << std::endl;
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

    // For all the 3d points in the hash_map
    for(const auto lm : sfm_data.GetLandmarks()) {
        const openMVG::Vec3 &pt = lm.second.X;
        positions.emplace_back(pt[0], pt[1], pt[2]);
    }

    std::vector<Alembic::Util::uint64_t> ids(positions.size());
    iota(begin(ids), end(ids), 0);

    OPoints partsOut(topObj, "particleShape1");
    OPointsSchema &pSchema = partsOut.getSchema();
    OPointsSchema::Sample psamp( move(V3fArraySample (positions)), move(UInt64ArraySample ( ids )));
    pSchema.set( psamp );

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
        openMVG::Mat3 R, K;
        openMVG::Vec3 t;
        openMVG::KRt_From_P(P, &K, &R, &t);
        
        // POSE
        // Compensate translation with rotation
        openMVG::Vec3 center = R.transpose() * t;
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

        // Get correct camera orientation
        M44d scale;
        scale[0][0] = -1;
        scale[1][1] = -1;
        scale[2][2] = -1;

        M44d scale2;
        scale2[0][0] = -1;
        scale2[1][1] = 1;
        scale2[2][2] = 1;
        
        xformMatrix = scale2*xformMatrix*scale;

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
        const float focalLengthPix = K(0,0);
        const float dx = K(0,2);
        const float dy = K(1,2);
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
