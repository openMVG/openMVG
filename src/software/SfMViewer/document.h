
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DOCUMENT
#define DOCUMENT

#include "software/SfM/SfMPinholeCamera.hpp"
using namespace openMVG;

#include "software/SfMViewer/Ply.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <vector>
#include <map>
#include <string>

struct Document
{
  std::vector<float> _vec_points;
  std::map<size_t, std::vector<size_t> > _map_visibility; //Inth camera see the Inth 3D point
  std::map<size_t, PinholeCamera > _map_camera;
  std::vector<std::string> _vec_imageNames;
  std::map<size_t, std::pair<size_t,size_t> > _map_imageSize;

  std::string _sDirectory;


  static bool readCamera(const std::string & sCam, PinholeCamera & cam)
  {
    std::vector<double> val;

    if (stlplus::extension_part(sCam) == "bin")
    {
      std::ifstream in(sCam.c_str(),
        std::ios::in|std::ios::binary);
      if (!in.is_open())	{
        std::cerr << "Error: failed to open file '" << sCam << "' for reading" << std::endl;
        return false;
      }
      val.resize(12);
      in.read((char*)&val[0],(std::streamsize)12*sizeof(double));
      if (in.fail())  {
        val.clear();
      }
    }
    else
      return false;

    if (val.size() == 3*4) //P Matrix
    {
      Mat34 P;
      P << val[0], val[3], val[6], val[9],
        val[1], val[4], val[7], val[10],
        val[2], val[5], val[8], val[11];

      Mat3 R,K;
      Vec3 t;
      KRt_From_P(P, &K, &R, &t);
      cam = PinholeCamera(K, R, t);
      return true;
    }
    return false;
  }


  bool load(const std::string & spath)
  {
    //-- Check if the required file are present.
    _sDirectory = spath;
    std::string sDirectoryPly = stlplus::folder_append_separator(_sDirectory) + "clouds";
    if (stlplus::is_file(stlplus::create_filespec(sDirectoryPly,"calib","ply"))
      && stlplus::is_file(stlplus::create_filespec(_sDirectory,"views","txt")))
    {
      //-- Read the ply file:
      Ply ply;
      size_t num_vertices = 0;
      // Read PLY header
      if (ply.open(stlplus::create_filespec(sDirectoryPly,"calib","ply")))
      {
        // ...
        std::cout << "PLY file opened";

        // Iterate over elements
        for (Ply::ElementsIterator it = ply.elements_begin();
          it != ply.elements_end(); ++it)
        {
          const Ply::Element& element = *it;
          if (element.name() != "vertex")
          {
            if (!ply.skip(element))
            {
              std::cerr << "Cannot skip element \"" << element.name() << '"' << std::endl;
              ply.close();
              return false;
            }
            continue;
          }
          num_vertices = element.count();

          //Reserve memory
          _vec_points.reserve(3*num_vertices);

          const size_t & num_elements = element.count();

          for (size_t i = 0; i != num_elements; ++i)
          {
            float pos[3];
            unsigned char color[3];
            float weight;
            std::vector<size_t> visibility;

            ply.read_begin(element);
            for (Ply::PropertiesIterator it2 =
              element.properties_begin();
              it2 != element.properties_end(); ++it2)
            {
              const Ply::Property& property = *it2;
              if (property.name() == "x")
                ply.read(property, pos[0]);
              else if (property.name() == "y")
                ply.read(property, pos[1]);
              else if (property.name() == "z")
                ply.read(property, pos[2]);
              else if (property.name() == "red")
                ply.read(property, color[0]);
              else if (property.name() == "green")
                ply.read(property, color[1]);
              else if (property.name() == "blue")
                ply.read(property, color[2]);
              else if (property.name() == "weight" || property.name() == "confidence")
                ply.read(property, weight);
              else if (property.name() == "visibility")
              {
                size_t visibility_count;
                ply.read_count(property, visibility_count);
                visibility.reserve(visibility_count);
                while (visibility_count--)
                {
                  int visibility_value;
                  ply.read_value(property, visibility_value);
                  visibility.push_back(visibility_value);
                  _map_visibility[visibility_value].push_back(i); //Jnth camera see the Inth point
                }
              }
              else if (!ply.skip(property))
              {
                std::cerr << "Cannot skip property \"" << property.name() << '"' << std::endl;
                ply.close();
                return false;
              }
            }
            ply.read_end();

            _vec_points.push_back(pos[0]);
            _vec_points.push_back(pos[1]);
            _vec_points.push_back(pos[2]);

            /*std::cout << '\n'
              << pos[0] <<' ' << pos[1] << ' ' << pos[2] << ' ';
            using namespace std;
            std::copy(visibility.begin(), visibility.end(), ostream_iterator<size_t>(std::cout, " "));
            */
          }
        }
      }
      ply.close();
    }
    else
    {
      std::cerr << "Required file(s) is missing" << std::endl;
    }

    // Read cameras
    std::string sDirectoryCam = stlplus::folder_append_separator(_sDirectory) + "cameras";
    //std::vector<std::string> vec_cameraNames = stlplus::folder_wildcard(sDirectoryCam,
    //  "*.bin", false, true);

    size_t camIndex = 0;
    //Read views file
    {
      std::ifstream iFilein(stlplus::create_filespec(_sDirectory,"views","txt").c_str());
      if (iFilein.is_open())
      {
        std::string temp;
        getline(iFilein,temp); //directory name
        getline(iFilein,temp); //directory name
        size_t nbImages;
        iFilein>> nbImages;
        while(iFilein.good())
        {
          getline(iFilein,temp);
          if (!temp.empty())
          {
            std::stringstream sStream(temp);
            std::string sImageName, sCamName;
            size_t w,h;
            float znear, zfar;
            sStream >> sImageName >> w >> h >> sCamName >> znear >> zfar;
            // Read the corresponding camera
            PinholeCamera cam;
            if (!readCamera(stlplus::folder_append_separator(sDirectoryCam) + sCamName, cam))
              std::cerr << "Cannot read camera" << std::endl;
            _map_camera[camIndex] = cam;

            _vec_imageNames.push_back(sImageName);
            _map_imageSize[camIndex] = std::make_pair(w,h);
            camIndex++;
          }
          temp.clear();
        }
      }
      std::cout << "\n Loaded image names : " << std::endl;
      std::copy(_vec_imageNames.begin(), _vec_imageNames.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
    }
    return !_map_camera.empty();
  }
};

#endif //DOCUMENT
