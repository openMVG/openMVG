// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/lm.hpp"

#include "testing/testing.h"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

#include <fstream>
#include <iostream>
#include <string>

using namespace openMVG;
using namespace svg;
using namespace std;

// Implementation of the problem found here:
// digilander.libero.it/foxes/optimiz/Optimiz1.htm

// Eigen LM functor to compute a mono-dimensional exponential regression
// We are looking for the exponential function that best fit the given point x,y
// f(x, c1, c2) = exp(c1*x) + c2
struct lm_Refine_functor : Functor<double>
{
  lm_Refine_functor(int inputs, int values,
    const Vec & x, const Vec & y): Functor<double>(inputs,values),
      _x(x), _y(y){}

  // The residual operator compute errors for the given x parameter
  int operator()(const Vec &x, Vec &fvec) const{
    // x contain the putative optimizer values
    double c1 = x[0];
    double c2 = x[1];

    // Evaluate the function cost for each input value
    for (Mat::Index i = 0; i < _x.rows(); ++i){
      fvec[i] = _y[i] - (exp(c1 * _x[i]) + c2);
    }
    return 0;
  }

  const Vec & _x, & _y; // Store data reference for cost evaluation
};

TEST(LM, MimimaSearchViaLM) {

  Vec x(10);
  x << .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
  Vec y(10);
  y << 1.7788, 1.6065, 1.4723, 1.3678, 1.2865, 1.2231, 1.1737, 1.1353, 1.1053, 1.0820;
  // We are looking for the exponential function that best fit the given point x,y
  // f(x, c1, c2) = exp(c1*x) + c2
  lm_Refine_functor functor(
    2,        // cardinal of the parameters: 2 (c1,c2)
    x.rows(),  // cardinal of computed residual
    x, y // Data that allow to compute the residual
  );
  using NumDiffT = Eigen::NumericalDiff<lm_Refine_functor>;
  NumDiffT numDiff(functor);

  Eigen::LevenbergMarquardt<NumDiffT > lm(numDiff);
  lm.parameters.maxfev = 100;

  // The starting point to optimize (a sufficiently approximate value)
  Vec xlm(2);
  xlm << 0.0, 0.0;

  // Optimization by using LevenbergMarquardt routine
  int info = lm.minimize(xlm);
  // Get back optimized value
  std::cout << "info" << info << std::endl;

  Vec minima = xlm;
  Vec2 GT(-.4999, .9999);
  EXPECT_MATRIX_NEAR(GT, xlm, 1e-3);

  // Evaluation of the residual of the found solution
  Vec fvec(10);
  functor(xlm, fvec);
  DOUBLES_EQUAL(0.0, fvec.norm(), 1e-3);
  // We cannot expect more precision since input data are limited in precision

  // Export computed result to a SVG file
  {
    //Draw input data and found curve (magnification by 10 for visualization):
    double dFactor = 10.0;
    svgDrawer svgSurface(6*dFactor, 4*dFactor);
    // Draw found curve
    Vec xFound(30),yFound(30);
    int cpt=0;
    for (double i=0.0; i < 6.0; i+=.2, ++cpt)  {
      xFound[cpt] = i*dFactor;
      yFound[cpt] = (4 - (exp(xlm[0] * i) + xlm[1]))*dFactor;
    }

    svgSurface.drawPolyline(
      xFound.data(), xFound.data()+30,
      yFound.data(), yFound.data()+30,
      svgStyle().stroke("blue", 1.f));

    //Draw point in a second time to put them in the top layer
    for (Vec::Index i = 0; i < x.rows(); ++i)
      svgSurface.drawCircle(x[i]*dFactor, (4-y[i])*dFactor, 1,
        svgStyle().fill("red").noStroke());

    std::ostringstream osSvg;
    osSvg << "exponentialRegression_unit_test.svg";
    std::ofstream svgFile( osSvg.str().c_str());
    svgFile << svgSurface.closeSvgFile().str();
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
