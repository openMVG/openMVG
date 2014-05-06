
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>

#include <GLFW/glfw3.h>

#include "software/SfMViewer/document.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/image/image.hpp"

#include "third_party/cmdLine/cmdLine.h"

static int running = 1;
static Document m_doc;
static int current_cam = -1;

static float x_offset=0;
static float y_offset=0;
static float z_offset=0;
static float normalized_focal = 1;

struct GLWImage {
    int width, height;
    GLuint texture;
    double opacity;
    int camera;
 };

//static GLWImage m_cur_image;
static std::vector< GLWImage > m_image_vector;

/* close callback */
void window_close_callback(GLFWwindow* window)
{
    running = 0;
}

void load_textures()
{
	int nbCams = m_doc._vec_imageNames.size();
	m_image_vector.resize(nbCams);

	for ( int i_cam=0; i_cam<nbCams; ++i_cam) {
		std::string sImageName = stlplus::create_filespec( stlplus::folder_append_separator(m_doc._sDirectory)+"images",
			m_doc._vec_imageNames[i_cam]);

		std::vector<unsigned char> img;
		int w,h,depth;
		if (ReadImage(sImageName.c_str(),	&img,	&w,	&h,	&depth)) {
			glEnable(GL_TEXTURE_2D);
			//std::cout << "Read image : " << sImageName << "\n" << std::endl;
			glDeleteTextures(1, &m_image_vector[i_cam].texture);

			// Create texture
			glGenTextures( 1, &m_image_vector[i_cam].texture);
			// select our current texture
		  glBindTexture(GL_TEXTURE_2D, m_image_vector[i_cam].texture);

		  m_image_vector[i_cam].width = w;
		  m_image_vector[i_cam].height = h;
		  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,  m_image_vector[i_cam].width,
		  		m_image_vector[i_cam].height, 0, GL_RGB, GL_UNSIGNED_BYTE,
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
  if( action != GLFW_PRESS ) {
  	return;
  }

  switch (k) {
  case GLFW_KEY_ESCAPE:
    running = 0;
    break;
  case GLFW_KEY_LEFT:
  	current_cam--;
    if (current_cam<0)  {
    	current_cam = m_doc._map_camera.size()-1;
    }
    break;
  case GLFW_KEY_RIGHT:
  	current_cam++;
    if (current_cam >= m_doc._map_camera.size())  {
    	current_cam = 0;
    }
    break;
  case GLFW_KEY_R:
  	x_offset=0;
  	y_offset= 0;
  	z_offset= 0;
  	normalized_focal = 1;
  	break;
  case GLFW_KEY_Q:
  	z_offset-= 0.1;
  	break;
  case GLFW_KEY_E:
  	z_offset+= 0.1;
    break;
  case GLFW_KEY_W:
  	y_offset+= 0.1;
  	break;
  case GLFW_KEY_S:
  	y_offset-= 0.1;
  	break;
  case GLFW_KEY_A:
  	x_offset+= 0.1;
  	break;
  case GLFW_KEY_D:
  	x_offset-= 0.1;
  	break;
  case GLFW_KEY_Z:
  	normalized_focal-= 0.1;
  	break;
  case GLFW_KEY_X:
  	normalized_focal+= 0.1;
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
    // convert opengl coordinates into the document information coordinates
    glPushMatrix();
    glMultMatrixf((GLfloat*)m_convert);

    // apply view offset
    openMVG::Mat4 offset_w = l2w_Camera(Mat3::Identity(), Vec3(x_offset,y_offset,z_offset));
    glMultMatrixd((GLdouble*)offset_w.data());

    // then apply current camera transformation
    const PinholeCamera & camera = m_doc._map_camera.find(current_cam)->second;
    openMVG::Mat4 l2w = l2w_Camera(camera._R, camera._t);

    glPushMatrix();
    glMultMatrixd((GLdouble*)l2w.data());

    glPointSize(3);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    //Draw Structure in GREEN (as seen from the current camera)
    size_t nbPoint = m_doc._vec_points.size()/3;
    size_t cpt = 0;
    glBegin(GL_POINTS);
    glColor3f(0.f,1.f,0.f);
    for(size_t i = 0; i < nbPoint; i++,cpt+=3) {
      glVertex3f(m_doc._vec_points[cpt], m_doc._vec_points[cpt+1], m_doc._vec_points[cpt+2]);
    }
    glEnd();

    glDisable(GL_CULL_FACE);

    for (int i_cam=0; i_cam<m_doc._vec_imageNames.size(); ++i_cam)
    {
    		const PinholeCamera & camera_i = m_doc._map_camera.find(i_cam)->second;

    		// Move frame to draw the camera i_cam by appliyin its inverse transformation
    		// Warning: translation has to be "fixed" to remove the current camera rotation

    		// Fix camera_i translation with current camera rotation inverse
    		Vec3 trans = camera._R.transpose() * camera_i._t;

    		// compute inverse transformation matrix from local to world
    		openMVG::Mat4 l2w_i = l2w_Camera(camera_i._R.transpose(), -trans);

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
   	    int w = m_doc._map_imageSize.find(i_cam)->second.first;
   	    int h = m_doc._map_imageSize.find(i_cam)->second.second;

   	    double focal = camera_i._K(0,0);
   	    // use principal point to adjust image center
   	    Vec2 pp(camera._K(0,2) , camera._K(1,2));

   	    Vec3 c1(    -pp[0]/focal * normalized_focal, (-pp[1]+h)/focal * normalized_focal, normalized_focal);
   	    Vec3 c2((-pp[0]+w)/focal * normalized_focal, (-pp[1]+h)/focal * normalized_focal, normalized_focal);
   	    Vec3 c3((-pp[0]+w)/focal * normalized_focal,     -pp[1]/focal * normalized_focal, normalized_focal);
   	    Vec3 c4(    -pp[0]/focal * normalized_focal,     -pp[1]/focal * normalized_focal, normalized_focal);

   	    // 2. Draw imagette
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

    glPopMatrix(); // go back to (document +offset) frame
    glPopMatrix(); // go back to identity
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
  	current_cam = 0;
    std::cout << "Press left or right key to navigate cameras ;-)" << std::endl;
    std::cout << "Move viewpoint with Q,W,E,A,S,D" << std::endl;
    std::cout << "Change Normalized focal with Z and X" << std::endl;
    std::cout << "Reset viewpoint position with R" << std::endl;
    std::cout << "Esc to quit" << std::endl;

  }
  else {
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

  window = glfwCreateWindow( 1000, 600, "SfmViewer", NULL, NULL );
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
