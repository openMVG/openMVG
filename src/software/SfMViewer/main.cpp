// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <iterator>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#ifdef _WIN32
#include <windows.h>
#endif
#endif

#include <GLFW/glfw3.h>

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image_io.hpp"
#include "third_party/progress/progress.hpp"

#include "third_party/cmdLine/cmdLine.h"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

static int running = 1;
static SfM_Data sfm_data;
static int current_cam = -1;

static float x_offset = 0.f;
static float y_offset = 0.f;
static float z_offset = 0.f;
static float normalized_focal = 1.f;

// Contiguous array of the valid camera index
static std::vector<IndexT> vec_cameras;

struct GLWImage {
    int width, height;
    GLuint texture;
 };

static std::vector< GLWImage > m_image_vector;

/* close callback */
void window_close_callback(GLFWwindow* window)
{
    running = 0;
}

void load_textures()
{
  const size_t nbCams = vec_cameras.size();
  m_image_vector.resize(nbCams);

  C_Progress_display my_progress_bar( nbCams, std::cout, "Textures loading, Please wait...\n" );
  for ( size_t i_cam=0; i_cam < nbCams; ++i_cam, ++my_progress_bar) {
    const View * view = sfm_data.GetViews().at(vec_cameras[i_cam]).get();
    const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);

  std::vector<unsigned char> img;
  int w,h,depth;
  if (ReadImage(srcImage.c_str(), &img, &w, &h, &depth)) {
    glEnable(GL_TEXTURE_2D);
    //std::cout << "Read image : " << sImageName << "\n" << std::endl;
    glDeleteTextures(1, &m_image_vector[i_cam].texture);

    // Create texture
    glGenTextures( 1, &m_image_vector[i_cam].texture);
    // select our current texture
    glBindTexture(GL_TEXTURE_2D, m_image_vector[i_cam].texture);

    m_image_vector[i_cam].width = w;
    m_image_vector[i_cam].height = h;
    glTexImage2D(GL_TEXTURE_2D, 0, (depth == 1) ? GL_LUMINANCE : GL_RGB,  m_image_vector[i_cam].width,
      m_image_vector[i_cam].height, 0, (depth == 1) ? GL_LUMINANCE : GL_RGB, GL_UNSIGNED_BYTE,
        &img[0]);

    glBindTexture(GL_TEXTURE_2D, m_image_vector[i_cam].texture);
  }
  }
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

void key(GLFWwindow* window, int k, int scancode, int action, int mod)
{
  if (action != GLFW_PRESS ) {
    return;
  }

  switch (k) {
  case GLFW_KEY_ESCAPE:
    running = 0;
    break;
  case GLFW_KEY_LEFT:
    --current_cam;
    if (current_cam < 0)  {
      current_cam = vec_cameras.size()-1;
    }
    break;
  case GLFW_KEY_RIGHT:
    ++current_cam;
    if (current_cam >= static_cast<int>(vec_cameras.size()))  {
      current_cam = 0;
    }
    break;
  case GLFW_KEY_R:
    x_offset = 0.f;
    y_offset = 0.f;
    z_offset = 0.f;
    normalized_focal = 1.f;
    break;
  case GLFW_KEY_Q:
    z_offset -= 0.1f;
    break;
  case GLFW_KEY_E:
    z_offset += 0.1f;
    break;
  case GLFW_KEY_W:
    y_offset += 0.1f;
    break;
  case GLFW_KEY_S:
    y_offset -= 0.1f;
    break;
  case GLFW_KEY_A:
    x_offset += 0.1f;
    break;
  case GLFW_KEY_D:
    x_offset -= 0.1f;
    break;
  case GLFW_KEY_KP_SUBTRACT:
    normalized_focal -= 0.1f;
    break;
  case GLFW_KEY_KP_ADD:
    normalized_focal += 0.1f;
    break;
  default:
    return;
  }
}

// the conversion matrix from OpenGL default coordinate system
//  to the camera coordinate system:
// [ 1  0  0  0] * [ x ] = [ x ]
//   0 -1  0  0      y      -y
//   0  0 -1  0      z      -z
//   0  0  0  1      1       1
const GLfloat m_convert[4][4] = {
  {1.f,  0.f,  0.f, 0.f},
  {0.f, -1.f,  0.f, 0.f},
  {0.f,  0.f, -1.f, 0.f},
  {0.f,  0.f,  0.f, 1.f}};

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
    // convert opengl coordinates into the document information coordinates
    glPushMatrix();
    glMultMatrixf((GLfloat*)m_convert);

    // apply view offset
    openMVG::Mat4 offset_w = l2w_Camera(Mat3::Identity(), Vec3(x_offset,y_offset,z_offset));
    glMultMatrixd((GLdouble*)offset_w.data());

    // then apply current camera transformation
    const View * view = sfm_data.GetViews().at(vec_cameras[current_cam]).get();
    const Pose3 pose = sfm_data.GetPoseOrDie(view);
    const openMVG::Mat4 l2w = l2w_Camera(pose.rotation(), pose.translation());

    glPushMatrix();
    glMultMatrixd((GLdouble*)l2w.data());

    glPointSize(3);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    //Draw Structure in GREEN (as seen from the current camera)
    glBegin(GL_POINTS);
    glColor3f(0.f,1.f,0.f);
    for (const auto & landmark_iter : sfm_data.GetLandmarks())
    {
      const Landmark & landmark = landmark_iter.second;
      glVertex3d(landmark.X(0), landmark.X(1), landmark.X(2));
    }
    glEnd();

    glDisable(GL_CULL_FACE);

    for (int i_cam=0; i_cam < static_cast<int>(vec_cameras.size()); ++i_cam)
    {
      const View * view = sfm_data.GetViews().at(vec_cameras[i_cam]).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      if (isPinhole(cam->getType()))
      {
        const Pinhole_Intrinsic * camPinhole = dynamic_cast<const Pinhole_Intrinsic*>(cam);

        // Move frame to draw the camera i_cam by applying its inverse transformation
        // Warning: translation has to be "fixed" to remove the current camera rotation

        // Fix camera_i translation with current camera rotation inverse
        const Vec3 trans = pose.rotation().transpose() * pose.translation();

        // compute inverse transformation matrix from local to world
        const openMVG::Mat4 l2w_i = l2w_Camera(pose.rotation().transpose(), -trans);
        // stack it and use it
        glPushMatrix();
        glMultMatrixd((GLdouble*)l2w_i.data());

        // 1. Draw optical center (RED) and image center (BLUE)
        glPointSize(3);
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_LIGHTING);

        glBegin(GL_POINTS);
        glColor3f(1.f,0.f,0.f);
        glVertex3f(0, 0, 0); // optical center
        glColor3f(0.f,0.f,1.f);
        glVertex3f(0, 0, normalized_focal); // image center
        glEnd();

        // compute image corners coordinated with normalized focal (f=normalized_focal)
        const int w = camPinhole->w();
        const int h = camPinhole->h();

        const double focal = camPinhole->focal();
        // use principal point to adjust image center
        const Vec2 pp = camPinhole->principal_point();

        const Vec3 c1(    -pp[0]/focal * normalized_focal, (-pp[1]+h)/focal * normalized_focal, normalized_focal);
        const Vec3 c2((-pp[0]+w)/focal * normalized_focal, (-pp[1]+h)/focal * normalized_focal, normalized_focal);
        const Vec3 c3((-pp[0]+w)/focal * normalized_focal,     -pp[1]/focal * normalized_focal, normalized_focal);
        const Vec3 c4(    -pp[0]/focal * normalized_focal,     -pp[1]/focal * normalized_focal, normalized_focal);

        // 2. Draw thumbnail
        if (i_cam == current_cam)
        {
          glEnable(GL_TEXTURE_2D);
          glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

          glBindTexture(GL_TEXTURE_2D, m_image_vector[i_cam].texture);

          glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
          glEnable(GL_BLEND);
          glDisable(GL_DEPTH_TEST);

          if (i_cam == current_cam) {
            glColor4f(0.5f,0.5f,0.5f, 0.7f);
          } else {
            glColor4f(0.5f,0.5f,0.5f, 0.5f);
          }

          glBegin(GL_QUADS);
          glTexCoord2d(0.0,1.0);    glVertex3d(c1[0], c1[1], c1[2]);
          glTexCoord2d(1.0,1.0);    glVertex3d(c2[0], c2[1], c2[2]);
          glTexCoord2d(1.0,0.0);    glVertex3d(c3[0], c3[1], c3[2]);
          glTexCoord2d(0.0,0.0);    glVertex3d(c4[0], c4[1], c4[2]);
          glEnd();

          glDisable(GL_TEXTURE_2D);
          glDisable(GL_BLEND); glEnable(GL_DEPTH_TEST);
        }

       // 3. Draw camera cone
        if (i_cam == current_cam) {
          glColor3f(1.f,1.f,0.f);
        } else {
          glColor3f(1.f,0.f,0.f);
        }
        glBegin(GL_LINES);
        glVertex3d(0.0,0.0,0.0); glVertex3d(c1[0], c1[1], c1[2]);
        glVertex3d(0.0,0.0,0.0); glVertex3d(c2[0], c2[1], c2[2]);
        glVertex3d(0.0,0.0,0.0); glVertex3d(c3[0], c3[1], c3[2]);
        glVertex3d(0.0,0.0,0.0); glVertex3d(c4[0], c4[1], c4[2]);
        glVertex3d(c1[0], c1[1], c1[2]); glVertex3d(c2[0], c2[1], c2[2]);
        glVertex3d(c2[0], c2[1], c2[2]); glVertex3d(c3[0], c3[1], c3[2]);
        glVertex3d(c3[0], c3[1], c3[2]); glVertex3d(c4[0], c4[1], c4[2]);
        glVertex3d(c4[0], c4[1], c4[2]); glVertex3d(c1[0], c1[1], c1[2]);
        glEnd();

        glPopMatrix(); // go back to current camera frame
      }
    }
    glPopMatrix(); // go back to (document +offset) frame
    glPopMatrix(); // go back to identity
  }
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename;
  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--sfmdata filename, the SfM_Data file to read]\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Read the SfM scene
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // List valid camera (view that have a pose & a valid intrinsic data)
  for(Views::const_iterator iter = sfm_data.GetViews().begin();
    iter != sfm_data.GetViews().end(); ++iter)
  {
    const View * view = iter->second.get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
      continue;

    vec_cameras.push_back(iter->first);
  }

  current_cam = 0;
  std::cout << "Press left or right key to navigate between cameras." << std::endl
    << "Move viewpoint with Q,W,E,A,S,D" << std::endl
    << "Change Normalized focal (camera cones size) with '+' and '-'" << std::endl
    << "Reset viewpoint position with R" << std::endl
    << "Esc to quit" << std::endl;

  //-- Create the GL window context
  GLFWwindow* window;
  int width, height;

  if ( !glfwInit() )
  {
    fprintf( stderr, "Failed to initialize GLFW\n" );
    exit( EXIT_FAILURE );
  }

  glfwWindowHint(GLFW_DEPTH_BITS, 16);

  window = glfwCreateWindow( 1000, 600, "SfmViewer", nullptr, nullptr );
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

  load_textures();

  // Main loop
  while (running)
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
