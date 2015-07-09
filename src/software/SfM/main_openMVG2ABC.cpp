#include <cstdlib>

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


int main(int argc, const char **argv) {
   
    //
    SfM_Data sfm_data;
    std::string filename = {"sfm_data.json"};// FIXME: give a correct file name
    if (!Load(sfm_data, filename, ESfM_Data(ALL))) {
        std::cerr << std::endl
        << "Error: The input project file \""<< filename << "\" cannot be read." << std::endl;
        return EXIT_FAILURE;
    }
  
    if (sfm_data.GetLandmarks().size() == 0) {
        std::cerr << "Error: no landmarks found" << std::endl;
        return EXIT_FAILURE;
    }

    // Fill vector with the values taken from OpenMVG 
    std::vector<V3f> positions;
    positions.reserve(sfm_data.GetLandmarks().size());

    // For all the 3d points in the hash_map
    for(const auto lm : sfm_data.GetLandmarks()) {
        const openMVG::Vec3 &pt = lm.second.X;
        positions.emplace_back(pt[0], pt[1], pt[2]);
    }

    std::vector<Alembic::Util::uint64_t> ids(positions.size());
    iota(begin(ids), end(ids), 0);

    // TODO: decide wether we want to keep HDF5 or switch to Ogawa 
    // FIXME : dynamic output filename instead of test
    OArchive archive(Alembic::AbcCoreHDF5::WriteArchive(), "test.abc");

    OObject topObj(archive, kTop);
    OPoints partsOut(topObj, "particleShape1");
    OPointsSchema &pSchema = partsOut.getSchema();

    OPointsSchema::Sample psamp( move(V3fArraySample (positions)), move(UInt64ArraySample ( ids )));
    pSchema.set( psamp );
    //=================================================================================
    //
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

        // Set the matrix in the sample
        XformSample xformsample;
        xformsample.setMatrix(xformMatrix);

        stringstream ss;
        ss << stlplus::basename_part(view->s_Img_path); 
        Alembic::AbcGeom::OXform xform(topObj, "camxform_" + ss.str());
        xform.getSchema().set(xformsample);

        // Camera intrinsic parameters
        // K is of the form:
        // 10861.9      -3.41061e-13 2879
        // 5.68434e-14  10861.9      1927.38
        // 0            0            1
        // Unit is pixels !!!
        // as we don't have the original film back, we set the horizontal aperture to 1mm
        // is equal to one and deduce the focal length from it
        const double widthPixels = static_cast<double>(cam->w());
        const double heightPixels = static_cast<double>(cam->h()); 
        const double focalLength = K(0,0)/widthPixels; // in mm
        const double horizontalAperture = 0.1; // in cm
        const double verticalAperture = 0.1*heightPixels/widthPixels; // in cm
        //const double verticalFilmOffset = 0.0;
        //const double horizontalFilmOffset = 0.0;
        OCamera camObj(xform, "camera_" + ss.str());
        CameraSample camSample;
        
        camSample.setFocalLength(focalLength);
        camSample.setHorizontalAperture(horizontalAperture);
        camSample.setVerticalAperture(verticalAperture);

        camObj.getSchema().set( camSample );
    } 

    return EXIT_SUCCESS;
}
