
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>

#include <GL/glfw3.h>

#include "software/SfMViewer/document.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/image/image.hpp"

#include "third_party/cmdLine/cmdLine.h"

static int running = 1;
static Document m_doc;
static int cur_cam = -1;

struct GLWImage {
    int width, height;
    GLuint texture;
    double opacity;
    int camera;
 };

static GLWImage m_cur_image;

/* close callback */
static int window_close_callback(GLFWwindow* window)
{
    running = 0;
    return GL_TRUE;
}

/* new window size */
void reshape( GLFWwindow* window, int width, int height )
{
  glViewport( 0, 0, (GLint) width, (GLint) height );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  GLfloat zNear = 1e-2;
  GLfloat zFar = 1e5;
  GLfloat aspect = float(width)/float(height);
  GLfloat fH = std::tan( float(60 / 360.0f * 3.14159f) ) * zNear;
  GLfloat fW = fH * aspect;
  glFrustum( -fW, fW, -fH, fH, zNear, zFar );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
}

void key( GLFWwindow* window, int k, int action)
{
  if( action != GLFW_PRESS ) return;

  bool bTextureChange = false;

  switch (k) {
  case GLFW_KEY_ESCAPE:
    running = 0;
    break;
  case GLFW_KEY_LEFT:
    cur_cam--;
    if (cur_cam<0)  {
      cur_cam = m_doc._map_camera.size()-1;
    }
    bTextureChange = true;
    break;
  case GLFW_KEY_RIGHT:
    cur_cam++;
    if (cur_cam >= m_doc._map_camera.size())  {
      cur_cam = 0;
    }
    bTextureChange = true;
    break;
  default:
    return;
  }

  if (bTextureChange)
  {
    std::string sImageName =
      stlplus::create_filespec(
        stlplus::folder_append_separator(m_doc._sDirectory)+"images",
        m_doc._vec_imageNames[cur_cam]);
    std::vector<unsigned char> img;
    int w,h,depth;
    if (ReadImage(sImageName.c_str(),
              &img,
              &w,
              &h,
              &depth))
    {
      glEnable(GL_TEXTURE_2D);
      std::cout << "Read image : " << sImageName << "\n" << std::endl;
      glDeleteTextures(1, &m_cur_image.texture);

      // Create texture
      glGenTextures( 1, &m_cur_image.texture);
      // select our current texture
      glBindTexture(GL_TEXTURE_2D, m_cur_image.texture);

      m_cur_image.width = w;
      m_cur_image.height = h;
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,  m_cur_image.width,
                m_cur_image.height, 0, GL_RGB, GL_UNSIGNED_BYTE,
                &img[0]);

      glBindTexture(GL_TEXTURE_2D, m_cur_image.texture);
    }
  }
}

// the conversion matrix from OpenGL default coordinate system
//  to the camera coordinate system:
// [ 1  0  0  0] * [ x ] = [ x ]
//   0 -1  0  0      y      -y
//   0  0 -1  0      z      -z
//   0  0  0  1      1       1
const GLfloat m_convert[4][4] = {
  {1.,  0.,  0., 0.},
  {0., -1.,  0., 0.},
  {0.,  0., -1., 0.},
  {0.,  0.,  0., 1.}};

//Local to World
static openMVG::Mat4 l2w_Camera(const Mat3 & R, const Vec3 & t)
{
  //World to Local
  /// given rotation matrix R and translation vector t,
  /// column-major matrix m is equal to:
  /// [ R11 R12 R13 t.x ]
  /// | R21 R22 R23 t.y |
  /// | R31 R32 R33 t.z |
  /// [ 0.0 0.0 0.0 1.0 ]

  //Local to World => Coordinates of the camera in the 3d space
  openMVG::Mat4 l2wmat = Mat4::Identity();
  l2wmat.block(0,0,3,4) = HStack(R,t);
  return l2wmat;
}

/* OpenGL draw function & timing */
static void draw(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  {
    const PinholeCamera & camera = m_doc._map_camera.find(cur_cam)->second;

    glPushMatrix();
    openMVG::Mat4 l2w = l2w_Camera(camera._R, camera._t);
    glMultMatrixf((GLfloat*)m_convert); // second, convert the camera coordinates to the opengl camera coordinates
    glMultMatrixd((GLdouble*)l2w.data());

    glPointSize(3);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    //glCallList(m_pointcloud);

    //Draw Structure
    size_t nbPoint = m_doc._vec_points.size()/3;
    size_t cpt = 0;
    glBegin(GL_POINTS);
    for(size_t i = 0; i < nbPoint; i++,cpt+=3) {
      glColor3f(0.f,1.f,0.f);
      glVertex3f(m_doc._vec_points[cpt], m_doc._vec_points[cpt+1], m_doc._vec_points[cpt+2]);
    }
    glEnd();


/*
    //-- Draw other cameras:
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND); glDisable(GL_DEPTH_TEST);

    glDisable(GL_CULL_FACE);
    for (std::map<size_t, SimpleCamera >::const_iterator iterCam = m_doc->map_camera.begin();
      iterCam != m_doc->map_camera.end(); ++iterCam)
    {
      glPushMatrix();

      if (std::distance(m_doc->map_camera.begin(),iterCam) != cur_cam)
      {
        Mat4 l2w = iterCam->second.l2w();
        //glMultMatrixd(l2w.data());
        glMultMatrixf((GLfloat*)m_convert); // second, convert the camera coordinates to the opengl camera coordinates
        //glMultMatrixf((GLfloat*)prj);       // first, project the points in the world coordinates to the camera coorrdinates
        glMultMatrixd((GLdouble*)l2w.data());
        int w = m_doc->map_imageSize.find(std::distance(m_doc->map_camera.begin(),iterCam))->second.first;
        int h = m_doc->map_imageSize.find(std::distance(m_doc->map_camera.begin(),iterCam))->second.second;
        double focal = iterCam->second.K(0,0);
        double maxx = 0.5*w/focal;
        double maxy = 0.5*h/focal;

        glBegin(GL_QUADS);
        glColor4f(0.5f,0.5f,0.5f,0.6f);
        glTexCoord2d(0.0,0.0);        glVertex3d(-maxx,-maxy,-1.0);
        glTexCoord2d(1.0,0.0);        glVertex3d(+maxx,-maxy,-1.0);
        glTexCoord2d(1.0,1.0);        glVertex3d(+maxx,+maxy,-1.0);
        glTexCoord2d(0.0,1.0);        glVertex3d(-maxx,+maxy,-1.0);
        glEnd();

      }
      glPopMatrix();
    }

    glDisable(GL_BLEND); glEnable(GL_DEPTH_TEST);
*/

    glPopMatrix();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND); glDisable(GL_DEPTH_TEST);

    //Draw the image

    int w = m_doc._map_imageSize.find(cur_cam)->second.first;
    int h = m_doc._map_imageSize.find(cur_cam)->second.second;
    double focal = camera._K(0,0);

    double maxx = 0.5 * w / focal;
    double maxy = 0.5 * h / focal;
    glEnable(GL_TEXTURE_2D);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, m_cur_image.texture);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.5f,0.5f,0.5f,0.6f);
    glBegin(GL_QUADS);
    glTexCoord2d(0.0,1.0);    glVertex3d(-maxx,-maxy,-1.0);
    glTexCoord2d(1.0,1.0);    glVertex3d(+maxx,-maxy,-1.0);
    glTexCoord2d(1.0,0.0);    glVertex3d(+maxx,+maxy,-1.0);
    glTexCoord2d(0.0,0.0);    glVertex3d(-maxx,+maxy,-1.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND); glEnable(GL_DEPTH_TEST);
  }
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Dir;
  cmd.add( make_option('i', sSfM_Dir, "sfmdir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << ' '
    << "[-i|--sfmdir SfM_Output path] "
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (m_doc.load(sSfM_Dir))
  {
    cur_cam = 0;
    std::cout << "Press left or right key to navigate ;-)" << std::endl;
    std::cout << "Esc to quit" << std::endl;
  }
  else{
    exit( EXIT_FAILURE);
  }

  //-- Create the GL window context
  GLFWwindow* window;
  int width, height;

  if( !glfwInit() )
  {
    fprintf( stderr, "Failed to initialize GLFW\n" );
    exit( EXIT_FAILURE );
  }

  glfwWindowHint(GLFW_DEPTH_BITS, 16);

  window = glfwCreateWindow( 600, 300, "SfmViewer", NULL, NULL );
  if (!window)
  {
    fprintf( stderr, "Failed to open GLFW window\n" );
    glfwTerminate();
    exit( EXIT_FAILURE );
  }

  // Set callback functions
  glfwSetWindowCloseCallback(window, window_close_callback);
  glfwSetWindowSizeCallback(window, reshape);
  glfwSetKeyCallback(window, key);

  glfwMakeContextCurrent(window);
  glfwSwapInterval( 1 );

  glfwGetWindowSize(window, &width, &height);
  reshape(window, width, height);

  m_cur_image.camera = -1;
  // Main loop
  while( running )
  {
    // Draw SfM Scene
    draw();

    // Swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  // Terminate GLFW
  glfwTerminate();

  // Exit program
  exit( EXIT_SUCCESS );
}
