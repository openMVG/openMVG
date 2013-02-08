 /*
 * @file main.cpp
 * @brief Example usage of svgDrawer.hpp functions
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

#include "svgDrawer.hpp"
using namespace svg;

int main(int argc, char * argv[])
{
  // Simple usage:
  {
    svgDrawer svgSurface; //Create a svg object
    // Add some drawing
    double S = 20.;
    for (double i = 0; i < 3.14*2; i+=.4) {
      double ax = cos(i)*S + S;
      double ay = sin(i)*S + S;
      double bx = cos(i+3.14/4.)*S + S;
      double by = sin(i+3.14/4.)*S + S;
      // Drawing are controlled by function
      //  svgStyle that configure the drawing option
      svgSurface.drawLine(ax, ay, bx, by,
        svgStyle().stroke("blue", 1));
    }
    //Export the SVG stream to a file
    std::string sFileName = "FirstExample.svg";
    std::ofstream svgFile(sFileName.c_str());
    svgFile << svgSurface.closeSvgFile().str();
    svgFile.close();
  }

  // Other drawing primitive:
  {
    svgDrawer svgSurface(20, 20); //Create a svg object
    // Add some drawing

    svgSurface.drawCircle(10,10, 4,
      svgStyle().stroke("red",1).fill("blue").tooltip("Hello"));

    svgSurface.drawSquare(4,4, 12, svgStyle().stroke("black"));

    svgSurface.drawText(8, 11, 6.f, "H", "green");

    //Export the SVG stream to a file
    std::string sFileName = "SecondExample.svg";
    std::ofstream svgFile(sFileName.c_str());
    svgFile << svgSurface.closeSvgFile().str();
    svgFile.close();
  }

  // Draw a cardiod with the svg polyline:
  // http://en.wikipedia.org/wiki/Cardioid
  {
    // Pre-compute (x,y) curve points
    size_t nbPoints = 120;
    std::vector<float> vec_x(nbPoints, 0.f), vec_y(nbPoints, 0.f);
    double S = 20.;
    for (size_t i = 0; i < nbPoints; ++i) {
      const double theta = i * 2 * M_PI / nbPoints;
      //-- Cardioid equation
      vec_x[i] = (3*S + S*(2.*sin(theta)-(sin(2.*theta))));
      vec_y[i] = (2*S - S*(2.*cos(theta)-(cos(2.*theta))));
    }
    // Create a svg surface and add the cardiod polyline
    svgDrawer svgSurface(6*S, 6*S); //Create a svg object
    svgSurface.drawPolyline(
      vec_x.begin(), vec_x.end(),
      vec_y.begin(), vec_y.end(),
      svgStyle().stroke("blue", 2));

    //Export the SVG stream to a file
    std::string sFileName = "ThirdExample.svg";
    std::ofstream svgFile(sFileName.c_str());
    svgFile << svgSurface.closeSvgFile().str();
    svgFile.close();
  }
  return EXIT_SUCCESS;
}
