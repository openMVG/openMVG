// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_DRAWING_HPP
#define OPENMVG_IMAGE_IMAGE_DRAWING_HPP

#include "openMVG/image/image_container.hpp"

namespace openMVG
{
namespace image
{

/**
* @brief Put the pixel in the image to the given color only if the point (xc,yc)
*  is inside the image.
* /!\ Be careful at the order (Y,X).
* @param yc Position on Y-axis
* @param xc Position on X-axis
* @param col Color to set
* @param[out] pim Output image where color is set
*/
template <typename Image, typename Color>
inline void SafePutPixel( int yc, int xc, const Color& col, Image *pim )
{
  if ( pim )
  {
    if ( pim->Contains( yc, xc ) )
    {
      ( *pim )( yc, xc ) = col;
    }
  }
}

/**
* @brief Bresenham approach to draw ellipse.
* @ref http://raphaello.univ-fcomte.fr/ig/algorithme/ExemplesGLUt/BresenhamEllipse.htm
* Add the rotation of the ellipse.
* As the algo. use symmetry we must use 4 rotations.
* @param xc Ellipse center on X-axis
* @param yx Ellipse center on Y-axis
* @param radiusA Main radius of ellipse
* @param radiusB Secondary radius of ellipse
* @param col Color of the ellipse
* @param[out] pim Output image
* @param angle Rotation angle of the Ellipse
*/
template <typename Image, typename Color>
void DrawEllipse( int xc, int yc, int radiusA, int radiusB,
                  const Color& col, Image *pim, double angle = 0.0 )
{
  int a = radiusA, b = radiusB;

  // Counter Clockwise rotation matrix.
  double matXY[4] = { cos( angle ), sin( angle ),
                      -sin( angle ), cos( angle )
                    };
  int x, y;
  double d1, d2;
  x = 0;
  y = b;
  d1 = b * b - a * a * b + a * a / 4;

  int rotX = ceil( matXY[0] * x + matXY[1] * y );
  int rotY = ceil( matXY[2] * x + matXY[3] * y );
  SafePutPixel( yc + rotY, xc + rotX, col, pim );
  rotX = matXY[0] * x - matXY[1] * y;
  rotY = matXY[2] * x - matXY[3] * y;
  SafePutPixel( yc + rotY, xc + rotX, col, pim );
  rotX = -matXY[0] * x - matXY[1] * y;
  rotY = -matXY[2] * x - matXY[3] * y;
  SafePutPixel( yc + rotY, xc + rotX, col, pim );
  rotX = -matXY[0] * x + matXY[1] * y;
  rotY = -matXY[2] * x + matXY[3] * y;
  SafePutPixel( yc + rotY, xc + rotX, col, pim );

  while ( a * a * ( y - .5 ) > b * b * ( x + 1 ) )
  {
    if ( d1 < 0 )
    {
      d1 += b * b * ( 2 * x + 3 );
      ++x;
    }
    else
    {
      d1 += b * b * ( 2 * x + 3 ) + a * a * ( -2 * y + 2 );
      ++x;
      --y;
    }
    rotX = matXY[0] * x + matXY[1] * y;
    rotY = matXY[2] * x + matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = matXY[0] * x - matXY[1] * y;
    rotY = matXY[2] * x - matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = -matXY[0] * x - matXY[1] * y;
    rotY = -matXY[2] * x - matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = -matXY[0] * x + matXY[1] * y;
    rotY = -matXY[2] * x + matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
  }
  d2 = b * b * ( x + .5 ) * ( x + .5 ) + a * a * ( y - 1 ) * ( y - 1 ) - a * a * b * b;
  while ( y > 0 )
  {
    if ( d2 < 0 )
    {
      d2 += b * b * ( 2 * x + 2 ) + a * a * ( -2 * y + 3 );
      --y;
      ++x;
    }
    else
    {
      d2 += a * a * ( -2 * y + 3 );
      --y;
    }
    rotX = matXY[0] * x + matXY[1] * y;
    rotY = matXY[2] * x + matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = matXY[0] * x - matXY[1] * y;
    rotY = matXY[2] * x - matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = -matXY[0] * x - matXY[1] * y;
    rotY = -matXY[2] * x - matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
    rotX = -matXY[0] * x + matXY[1] * y;
    rotY = -matXY[2] * x + matXY[3] * y;
    SafePutPixel( yc + rotY, xc + rotX, col, pim );
  }
}


/**
* @brief Bresenham approach do not allow to draw concentric circle without holes.
* So it's better the use the Andres method.
* @ref http://fr.wikipedia.org/wiki/Algorithme_de_tracé_de_cercle_d'Andres.
* @param x Circle center on X-axis
* @param y Circle center on Y-axis
* @param col Color of the circle
* @param[out] pim Output image
*/
template <typename Image, typename Color>
void DrawCircle( int x, int y, int radius, const Color& col, Image *pim )
{
  Image &im = *pim;
  if (  im.Contains( y + radius, x + radius )
        || im.Contains( y + radius, x - radius )
        || im.Contains( y - radius, x + radius )
        || im.Contains( y - radius, x - radius ) )
  {
    int x1 = 0;
    int y1 = radius;
    int d = radius - 1;
    while ( y1 >= x1 )
    {
      // Draw the point for each octant.
      SafePutPixel( y1 + y,  x1 + x, col, pim );
      SafePutPixel( x1 + y,  y1 + x, col, pim );
      SafePutPixel( y1 + y, -x1 + x, col, pim );
      SafePutPixel( x1 + y, -y1 + x, col, pim );
      SafePutPixel( -y1 + y,  x1 + x, col, pim );
      SafePutPixel( -x1 + y,  y1 + x, col, pim );
      SafePutPixel( -y1 + y, -x1 + x, col, pim );
      SafePutPixel( -x1 + y, -y1 + x, col, pim );
      if ( d >= 2 * x1 )
      {
        d = d - 2 * x1 - 1;
        x1 += 1;
      }
      else
      {
        if ( d <= 2 * ( radius - y1 ) )
        {
          d = d + 2 * y1 - 1;
          y1 -= 1;
        }
        else
        {
          d = d + 2 * ( y1 - x1 - 1 );
          y1 -= 1;
          x1 += 1;
        }
      }
    }
  }
}


/**
* @brief Bresenham algorithm used to draw a line
* @param xa Start of the line on X-axis
* @param ya Start of the line on Y-axis
* @param xb End of the line on X-axis
* @param yb End of the line on Y-axis
* @param col Color of the line
* @param[out] pim Output image
*/
template <typename Image, typename Color>
void DrawLine( int xa, int ya, int xb, int yb, const Color& col, Image *pim )
{
  Image &im = *pim;

  // If one point is outside the image
  // Replace the outside point by the intersection of the line and
  // the limit (either x=width or y=height).
  if ( !im.Contains( ya, xa ) || !im.Contains( yb, xb ) )
  {

    int width = pim->Width(), height = pim->Height();
    const bool xdir = xa < xb, ydir = ya < yb;
    float nx0 = float( xa ), nx1 = float( xb ), ny0 = float( ya ), ny1 = float( yb ),
          &xleft = xdir ? nx0 : nx1,  &yleft = xdir ? ny0 : ny1,
           &xright = xdir ? nx1 : nx0, &yright = xdir ? ny1 : ny0,
            &xup = ydir ? nx0 : nx1,    &yup = ydir ? ny0 : ny1,
             &xdown = ydir ? nx1 : nx0,  &ydown = ydir ? ny1 : ny0;

    if ( xright < 0 || xleft >= width )
    {
      return;
    }
    if ( xleft < 0 )
    {
      yleft -= xleft * ( yright - yleft ) / ( xright - xleft );
      xleft  = 0;
    }
    if ( xright >= width )
    {
      yright -= ( xright - width ) * ( yright - yleft ) / ( xright - xleft );
      xright  = float( width ) - 1;
    }
    if ( ydown < 0 || yup >= height )
    {
      return;
    }
    if ( yup < 0 )
    {
      xup -= yup * ( xdown - xup ) / ( ydown - yup );
      yup  =  0;
    }
    if ( ydown >= height )
    {
      xdown -= ( ydown - height ) * ( xdown - xup ) / ( ydown - yup );
      ydown  =  float( height ) - 1;
    }

    xa = ( int ) xleft;
    xb = ( int ) xright;
    ya = ( int ) yleft;
    yb = ( int ) yright;
  }

  int xbas, xhaut, ybas, yhaut;
  // Check the condition ybas < yhaut.
  if ( ya <= yb )
  {
    xbas = xa;
    ybas = ya;
    xhaut = xb;
    yhaut = yb;
  }
  else
  {
    xbas = xb;
    ybas = yb;
    xhaut = xa;
    yhaut = ya;
  }
  // Initialize slope.
  int x, y, dx, dy, incrmX, incrmY, dp, N, S;
  dx = xhaut - xbas;
  dy = yhaut - ybas;
  if ( dx > 0 ) // If xhaut > xbas we will increment X.
  {
    incrmX = 1;
  }
  else
  {
    incrmX = -1; // else we will decrement X.
    dx *= -1;
  }
  if ( dy > 0 ) // Positive slope will increment X.
  {
    incrmY = 1;
  }
  else          // Negative slope.
  {
    incrmY = -1;
  }
  if ( dx >= dy )
  {
    dp = 2 * dy - dx;
    S = 2 * dy;
    N = 2 * ( dy - dx );
    y = ybas;
    x = xbas;
    while ( x != xhaut )
    {
      SafePutPixel( y, x, col, pim );
      x += incrmX;
      if ( dp <= 0 ) // Go in direction of the South Pixel.
      {
        dp += S;
      }
      else           // Go to the North.
      {
        dp += N;
        y += incrmY;
      }
    }
  }
  else
  {
    dp = 2 * dx - dy;
    S = 2 * dx;
    N = 2 * ( dx - dy );
    x = xbas;
    y = ybas;
    while ( y < yhaut )
    {
      SafePutPixel( y, x, col, pim );
      y += incrmY;
      if ( dp <= 0 ) // Go in direction of the South Pixel.
      {
        dp += S;
      }
      else           // Go to the North.
      {
        dp += N;
        x += incrmX;
      }
    }
  }
  SafePutPixel( y, x, col, pim );
}


/**
* @brief Draw filled circle
* Exterior point computed with bresenham approach
* @see DrawCircle
* @param x Circle center on X-axis
* @param y Circle center on Y-axis
* @param col Color of the circle
* @param[out] pim Output image
*/
template <typename Image, typename Color>
void FilledCircle( int x, int y, int radius, const Color& col, Image *pim )
{
  Image &im = *pim;
  if (  im.Contains( y + radius, x + radius )
        || im.Contains( y + radius, x - radius )
        || im.Contains( y - radius, x + radius )
        || im.Contains( y - radius, x - radius ) )
  {
    int x1 = 0;
    int y1 = radius;
    int d = radius - 1;
    while ( y1 >= x1 )
    {
      DrawLine( x1 + x, y1 + y, x1 + x, -y1 + y, col, pim );
      DrawLine( y1 + x, x1 + y, y1 + x, -x1 + y, col, pim );
      DrawLine( -x1 + x, y1 + y, -x1 + x, -y1 + y, col, pim );
      DrawLine( -y1 + x, x1 + y, -y1 + x, -x1 + y, col, pim );
      if ( d >= 2 * x1 )
      {
        d = d - 2 * x1 - 1;
        x1 += 1;
      }
      else
      {
        if ( d <= 2 * ( radius - y1 ) )
        {
          d = d + 2 * y1 - 1;
          y1 -= 1;
        }
        else
        {
          d = d + 2 * ( y1 - x1 - 1 );
          y1 -= 1;
          x1 += 1;
        }
      }
    }
  }
}


/**
* @brief Draw a serie of circles along the line, the algorithm is slow but accurate
* @param xa Start of the line on X-axis
* @param ya Start of the line on Y-axis
* @param xb End of the line on X-axis
* @param yb End of the line on Y-axis
* @param col Color of the line
* @param thickness Thiness of the line to draw
* @param[out] pim Output image
*/
template <typename Image, typename Color>
void DrawLineThickness( int xa, int ya, int xb, int yb, const Color& col, int thickness, Image *pim )
{
  Image &im = *pim;
  int halfThickness = ( thickness + 1 ) / 2;

  // If one point is outside the image
  // Replace the outside point by the intersection of the line and
  // the limit (either x=width or y=height).
  if ( !im.Contains( ya, xa ) || !im.Contains( yb, xb ) )
  {

    int width = pim->Width(), height = pim->Height();
    const bool xdir = xa < xb, ydir = ya < yb;
    float nx0 = float( xa ), nx1 = float( xb ), ny0 = float( ya ), ny1 = float( yb ),
          &xleft = xdir ? nx0 : nx1,  &yleft = xdir ? ny0 : ny1,
           &xright = xdir ? nx1 : nx0, &yright = xdir ? ny1 : ny0,
            &xup = ydir ? nx0 : nx1,    &yup = ydir ? ny0 : ny1,
             &xdown = ydir ? nx1 : nx0,  &ydown = ydir ? ny1 : ny0;

    if ( xright < 0 || xleft >= width )
    {
      return;
    }
    if ( xleft < 0 )
    {
      yleft -= xleft * ( yright - yleft ) / ( xright - xleft );
      xleft  = 0;
    }
    if ( xright >= width )
    {
      yright -= ( xright - width ) * ( yright - yleft ) / ( xright - xleft );
      xright  = float( width ) - 1;
    }
    if ( ydown < 0 || yup >= height )
    {
      return;
    }
    if ( yup < 0 )
    {
      xup -= yup * ( xdown - xup ) / ( ydown - yup );
      yup  =  0;
    }
    if ( ydown >= height )
    {
      xdown -= ( ydown - height ) * ( xdown - xup ) / ( ydown - yup );
      ydown  =  float( height ) - 1;
    }

    xa = ( int ) xleft;
    xb = ( int ) xright;
    ya = ( int ) yleft;
    yb = ( int ) yright;
  }

  int xbas, xhaut, ybas, yhaut;
  // Check the condition ybas < yhaut.
  if ( ya <= yb )
  {
    xbas = xa;
    ybas = ya;
    xhaut = xb;
    yhaut = yb;
  }
  else
  {
    xbas = xb;
    ybas = yb;
    xhaut = xa;
    yhaut = ya;
  }
  // Initialize slope.
  int x, y, dx, dy, incrmX, incrmY, dp, N, S;
  dx = xhaut - xbas;
  dy = yhaut - ybas;
  if ( dx > 0 ) // If xhaut > xbas we will increment X.
  {
    incrmX = 1;
  }
  else
  {
    incrmX = -1; // else we will decrement X.
    dx *= -1;
  }
  if ( dy > 0 ) // Positive slope will increment X.
  {
    incrmY = 1;
  }
  else          // Negative slope.
  {
    incrmY = -1;
  }
  if ( dx >= dy )
  {
    dp = 2 * dy - dx;
    S = 2 * dy;
    N = 2 * ( dy - dx );
    y = ybas;
    x = xbas;
    while ( x != xhaut )
    {
      DrawCircle( x, y, halfThickness, col, pim );
      x += incrmX;
      if ( dp <= 0 ) // Go in direction of the South Pixel.
      {
        dp += S;
      }
      else           // Go to the North.
      {
        dp += N;
        y += incrmY;
      }
    }
  }
  else
  {
    dp = 2 * dx - dy;
    S = 2 * dx;
    N = 2 * ( dx - dy );
    x = xbas;
    y = ybas;
    while ( y < yhaut )
    {
      DrawCircle( x, y, halfThickness, col, pim );
      y += incrmY;
      if ( dp <= 0 ) // Go in direction of the South Pixel.
      {
        dp += S;
      }
      else           // Go to the North.
      {
        dp += N;
        x += incrmX;
      }
    }
  }
  DrawCircle( x, y, halfThickness, col, pim );
}

} // namespace image
} // namespace openMVG

#endif  // OPENMVG_IMAGE_IMAGE_DRAWING_HPP
