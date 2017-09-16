// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef CGLWINDOW_HPP
#define CGLWINDOW_HPP

//include header file for glfw library so that we can use OpenGL
#if defined(__APPLE__)
#include <OpenGL/gl.h>
#endif
#include <GLFW/glfw3.h>

/// Basic class to handle a GLFW window (no interaction)
/// GLFW have be initialized before calling this class
struct CGlWindow
{
  bool _bRunning;

  GLFWwindow * _window;
  int _width, _height;

  CGlWindow() : _window(NULL), _width(-1), _height(-1) { }

  bool Init(int w, int h, const std::string & sWindowName)
  {
    glfwWindowHint(GLFW_DEPTH_BITS, 16);

    _width = w;
    _height = h;

    _window = glfwCreateWindow( w, h, sWindowName.c_str(), NULL, NULL );
    if (!_window)
    {
      glfwTerminate();
      return false;
    }
    else
    {
      _bRunning = true;
    }

    glfwMakeContextCurrent(_window);
    glfwSwapInterval( 1 );

    return true;
  }

  void SetOrtho(int W, int H)
  {
    //Draw full screen image
    glViewport( 0, 0, _width, _height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, W, H, 0.0, -1, 1);
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
  }

  void DrawFullScreenTexQuad(int W, int H)
  {
    // Draw a fullscreen Image Quad
    glBegin(GL_QUADS);
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glTexCoord2f(0.0f, 0.0f); glVertex3f(0.f, 0.f, 0.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex3f(GLfloat(W), 0.0f, 0.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(GLfloat(W), GLfloat(H), 0.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(0.f, GLfloat(H), 0.0f);
    glEnd();
  }

  void Swap()
  {
    glfwSwapBuffers(_window);
    glfwPollEvents(); // Poll for and process events
  }
};


#endif // CGLWINDOW_HPP
