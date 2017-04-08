// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
// Copyright (c) 2016 nomoko AG, Srivathsan Murali<srivathsan@nomoko.camera>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef _DOMSET_INTERFACE_H_
#define _DOMSET_INTERFACE_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Dense>
#include <types.h>

namespace nomoko{
  struct Interface {
    std::vector<Camera> cameras;
    void addCamera(const Eigen::Matrix3f& _K) {
      Camera c;
      c.K = _K;
      cameras.push_back(c);
    }

    std::vector<View> views;
    void addView(
        const Eigen::Matrix3f& _rot,
        const Eigen::Vector3f& _trans,
        const unsigned int& _cameraId,
        const std::string& _filename) {
      View v;
      v.rot = _rot;
      v.trans = _trans;
      v.cameraId = _cameraId;
      v.filename = _filename;
      views.push_back(v);
    }

    std::vector<Point> points;
    void addPoint(
        const Eigen::Vector3f& _pos,
        const std::vector<unsigned int> & _viewList) {
      Point p;
      p.pos = _pos;
      p.viewList = _viewList;
      points.push_back(p);
    }

    bool save(const std::string& domsetFilePath) {
      std::ofstream df (domsetFilePath);
      if(!df.is_open()){
        std::cout<< "Cant open file " << domsetFilePath << std::endl;
        return false;
      }

      df << "Domset\n";
      df << cameras.size() << " "
        << views.size() << " "
        << points.size() << std::endl;

      for(const Camera c : cameras) {
        df << c.K << std::endl;
        df << c.width << " " << c.height << std::endl;
      }

      for(const View v : views) {
        df << v.cameraId << std::endl
          << v.rot << std::endl
          << v.trans.transpose() << std::endl
          << v.filename <<std::endl;
      }

      for (const Point p : points) {
        df << p.pos.transpose() << std::endl;
        df << p.viewList.size();
        for(const unsigned int vId : p.viewList) {
          df << " " << vId;
        }
        df << std::endl;
      }
      df.close();
      return true;
    }

    bool load(const std::string& domsetFilepath) {
      std::ifstream df (domsetFilepath);
      if(!df.is_open()) {
        std::cout<< "Cant open file " << domsetFilepath << std::endl;
        exit(0);
      }

      std::string header;
      df >> header;

      int tmp;
      df >> tmp;
      const int numCameras = tmp;
      df >> tmp;
      const int numViews= tmp;
      df >> tmp;
      const int numPoints= tmp;

      cameras.reserve(numCameras);
      views.reserve(numViews);
      points.reserve(numPoints);

      for(int i =0; i < numCameras; i++ ) {
        Camera c;
        df >> c.K(0,0) >> c.K(0,1) >> c.K(0,2)
          >> c.K(1,0) >> c.K(1,1) >> c.K(1,2)
          >> c.K(2,0) >> c.K(2,1) >> c.K(2,2);
        df >> c.width >> c.height;
        cameras.push_back(c);
      }

      for(int i =0; i < numViews; i++ ) {
        View v;
        df >> v.cameraId;
        df >> v.rot(0,0) >> v.rot(0,1) >> v.rot(0,2)
          >> v.rot(1,0) >> v.rot(1,1) >> v.rot(1,2)
          >> v.rot(2,0) >> v.rot(2,1) >> v.rot(2,2);
        df >> v.trans(0) >> v.trans(1) >> v.trans(2);
        df >> v.filename;
        views.push_back(v);
      }
      int vSize;
      int tmpV;
      for(int i =0; i < numPoints; i++ ) {
        Point p;
        df >> p.pos(0) >> p.pos(1) >> p.pos(2);
        df >> vSize;
        for(int j=0; j < vSize; j++){
          df >> tmpV;
          p.viewList.push_back(tmpV);
        }
        points.push_back(p);
      }
    }

    void print() {
      std::stringstream ss;
      ss << "Domset\n";
      for(const Camera c : cameras) {
        ss << c.K << std::endl;
      }

      for(const View v : views) {
        ss << v.cameraId << std::endl
          << v.rot << std::endl
          << v.trans.transpose() << std::endl
          << v.filename <<std::endl;
      }

      for (const Point p : points) {
        ss<< p.pos.transpose() << std::endl;
        ss<< p.viewList.size();
        for(const unsigned int vId : p.viewList) {
          ss<< " " << vId;
        }
        ss<< std::endl;
      }
      std::cerr<< ss.str() << std::endl;
    }
  }; // struct DomsetInterface
} // namespace nomoko
#endif // _DOMSET_INTERFACE_H_
