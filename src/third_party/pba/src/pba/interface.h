////////////////////////////////////////////////////////////////////////////
//	File:		    interface.h
//	Author:		    Changchang Wu (ccwu@cs.washington.edu)
//	Description :   c-like interface of the library
//
//  Copyright (c) 2011  Changchang Wu (ccwu@cs.washington.edu)
//    and the University of Washington at Seattle 
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation; either
//  Version 3 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
using std::cerr; 
#include "pba.h"

////////////////////////////////
//This file give several template interfaces similar to other softwares
//      run_sfm_pba(...)                     an interface similar  to Snavely's "run_sfm"
//      run_pba(motstruct, vmask,...)        an interface similar to Lourakis's "sba"


#define CHECK_PARAM_SUPPORT(param) \
    if(param) {std::cerr << #param << " is not supported\n"; return false; }


template<class camera_params_t, class v3_t, class Float>
bool run_sfm_pba(int num_pts, int num_cameras, int ncons,
        char *vmask, Float *projections,
        int est_focal_length, int const_focal_length, //? duplicate parameters?
        int undistort,  int explicit_camera_centers,
        camera_params_t *init_camera_params, v3_t *init_pts,
        int use_constraints,   int use_point_constraints,
        v3_t *pt_constraints,  Float pt_constraint_weight,  int fix_points,
        int optimize_for_fisheye, Float eps2,  Float *Vout, Float *Sout,
        Float *Uout, Float *Wout, int max_iteration = 150)
{
    ///////////////////////////////////////////////
    if(vmask == NULL || projections == NULL)  return false;
    CHECK_PARAM_SUPPORT(explicit_camera_centers == false);
    CHECK_PARAM_SUPPORT(optimize_for_fisheye);
    CHECK_PARAM_SUPPORT(use_constraints);
    CHECK_PARAM_SUPPORT(use_point_constraints);
    CHECK_PARAM_SUPPORT(fix_points);
    /////////////////////

    vector<CameraT>        camera_data(num_cameras);    //camera (input/ouput)
    vector<Point3D>        point_data(num_pts);     //3D point(iput/output)
    vector<Point2D>        measurements;   //measurment/projection vector
    vector<int>            camidx, ptidx;  //index of camera/point for each projection

    for(int i = 0, k = 0; i < num_pts; ++i)
    {
        point_data[i].SetPoint(init_pts[i].p);
        for(int j = 0; j < num_cameras; ++j)
        {
            if(vmask[i * num_cameras + j] == 0) continue;
            measurements.push_back(Point2D(-projections[2 * k], -projections[2 * k + 1]));
            camidx.push_back(j);
            ptidx.push_back(i);
            ++k;
        }
    }

    for(int i = 0; i < num_cameras; ++i)
    {
        CameraT& cam = camera_data[i];
        cam.SetMatrixRotation(init_camera_params[i].R);
        cam.SetCameraCenterAfterRotation(init_camera_params[i].t);
        cam.SetFocalLength(init_camera_params[i].f);
        cam.SetProjectionDistortion(init_camera_params[i].k[0]);
    }

    std::cout << "[PBA] ncam = " << num_cameras << "; npt = " << num_pts << "; nproj = " << ptidx.size() << '\n';

    static ParallelBA gpba;
    gpba.GetInternalConfig()->__lm_max_iteration = max_iteration;
    gpba.GetInternalConfig()->__lm_damping_auto_switch = 2.0f;
    gpba.EnableRadialDistortion(undistort? ParallelBA::PBA_PROJECTION_DISTORTION: ParallelBA::PBA_NO_DISTORTION);
    gpba.SetCameraData(camera_data.size(),  &camera_data[0]);                        //set camera parameters
    gpba.SetPointData(point_data.size(), &point_data[0]);                            //set 3D point data
    gpba.SetProjection(measurements.size(), &measurements[0], &ptidx[0], &camidx[0]);//set the projections
    gpba.RunBundleAdjustment();

    /////////save the new  camera parameters
    for(int i = 0; i < num_cameras; ++i)
    {
        CameraT& cam = camera_data[i];
        cam.GetMatrixRotation(init_camera_params[i].R);
        cam.GetCameraCenter(init_camera_params[i].t);
        init_camera_params[i].f = cam.GetFocalLength();
        init_camera_params[i].k[0] = cam.GetProjectionDistortion();
        init_camera_params[i].k[1] = 0;
    }

    ////save new point data
    for(int i = 0; i < num_pts; ++i)
    {
        point_data[i].GetPoint(init_pts[i].p);
    }
    return true;
}




template <class Float> 
bool run_pba(int num_cameras, int num_pts, int numprojs, Float *motstruct, const Float*imgpts, const char*vmask )
{
    //suppose each camera has F(1), R(3x3), T(3)
    const int cam_param_num = 13;
    Float * cam = motstruct;
    Float * pts = motstruct + cam_param_num * num_cameras;

    /////////////////////////////////////////////
    vector<CameraT>        camera_data(num_cameras);    //camera (input/ouput)
    vector<Point3D>        point_data(num_pts);     //3D point(iput/output)
    vector<Point2D>        measurements(numprojs);   //measurment/projection vector
    vector<int>            camidx, ptidx;  //index of camera/point for each projection

    for(int i = 0; i < num_cameras; ++i, cam += cam_param_num)
    {
        CameraT& camera = camera_data[i];
        camera.SetFocalLength(cam[0]);
        camera.SetMatrixRotation(cam + 1);
        camera.SetTranslation(cam + 10);
    }

    for(int i = 0; i < num_pts; ++i, pts += 3 )
    {
        point_data[i].SetPoint(pts);
    }

    for(int i = 0; i < numprojs; ++i, imgpts += 2)
    {
        measurements[i].SetPoint2D(imgpts[0], imgpts[1]);
    }

    camidx.reserve(numprojs);
    ptidx.reserve(numprojs);
    for (int j = 0; j < num_pts; ++j)
    {
        for(int i = 0; i < num_cameras; ++i, ++vmask)
        {
            if(*vmask == 0) continue;
            camidx.push_back(i);
            ptidx.push_back(j);
        }
    }

    //////////////////////////////////////
    static ParallelBA gpba;
    gpba.SetCameraData(camera_data.size(),  &camera_data[0]);                        //set camera parameters
    gpba.SetPointData(point_data.size(), &point_data[0]);                            //set 3D point data
    gpba.SetProjection(measurements.size(), &measurements[0], &ptidx[0], &camidx[0]);//set the projections
    gpba.RunBundleAdjustment();


    //save new camera data
    cam = motstruct;
    for(int i = 0; i < num_cameras; ++i, cam += cam_param_num)
    {
        CameraT& camera = camera_data[i];
        cam[0] = camera.GetFocalLength();
        camera.GetMatrixRotation(cam + 1);
        camera.GetTranslation(cam + 10);
    }
 
    ////save new point data
    pts = motstruct + cam_param_num * num_cameras;
    for(int i = 0; i < num_pts; ++i, pts += 3 )
    {
        point_data[i].GetPoint(pts);
    }
    return true;
}