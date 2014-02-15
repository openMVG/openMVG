/*
 * This script export a full website using WebGL for visualization of openMVG output 
 * 
 * it only needs the user to provide input folder (ie SfM_Output folder of openMVG)
 * and an output folder (that will be created ou overwrited so be careful)
 *
 * If you want custom web site for your project please do as follow :
 * 1 - use this script to generate a standard website
 * 2 - customize anything you want (hmtl/css/js) except sceneData.js
 *    ( dont change variables names or directory structure)  
 * 3 - save your site directory structure in a location (ex A)
 * 4 - export a new openMVG project 
 * 5 - copy sceneData.js to your folder (ie A/scripts)
 * 6 - copy images folder content into your images folder (ie A/images)
 * 7 Enjoy ! 
 */


#include <iostream>
#include <cstdlib>
#include <iomanip>


#include "software/SfMViewer/document.h"
#include "third_party/cmdLine/cmdLine.h"
#include "openMVG/image/image.hpp"




/**
  * Prepare WebGL output folder 
  * {WebGL}/style              -- CSS Style file 
  * {WebGL}/scripts            -- js scripts
  * {WebGL}/images             -- pictures 
  * @param sOutputDir base folder where project will be saved  
*/

void prepareFolders( const std::string & sOutputDir )
{
  // Create initial folder 
  if( ! stlplus::folder_exists(sOutputDir) )
  {
    stlplus::folder_create( sOutputDir ) ;
  }

	// Create style folder
  if( ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir ) + "style" ) )
  {
    stlplus::folder_create( stlplus::folder_append_separator( sOutputDir ) + "style" ) ;
  }

	// Create scripts folder 
  if( ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir ) + "scripts" ) )
  {
    stlplus::folder_create( stlplus::folder_append_separator( sOutputDir ) + "scripts" ) ;
  }

	// Create image folder
  if( ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir ) + "images" ) )
  {
    stlplus::folder_create( stlplus::folder_append_separator( sOutputDir ) + "images" ) ;
  }


	// Check validity 
  if( ! stlplus::folder_exists( sOutputDir ) || 
    ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir) + "style" ) || 
    ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir) + "scripts" ) ||
    ! stlplus::folder_exists( stlplus::folder_append_separator( sOutputDir) + "images" ) )
  {
    std::cerr << "Cannot create project folder" << std::endl ; 
    exit( EXIT_FAILURE ) ; 
  }	
}


/**
  * Compute new width and height in order to respect original ratio 
  * @param newWidth (in: max new width) ( out: computed new width)
  * @param newHeight (in : max new height) ( out: computed new height)
  * @param oldWidth original image width 
  * @param oldHeight original image height 
  */
static inline void RessamplingSize( int & newWidth , int & newHeight , const int oldWidth , const int oldHeight )
{
	if( oldWidth > oldHeight )
	{
		const double ratio = (double) newWidth / (double) oldWidth ; 
		newHeight = ratio * (double) oldHeight ; 
	}
	else
	{
		const double ratio = (double) newHeight / (double) oldHeight ;
		newWidth = ratio * (double) oldWidth ; 
	}
}

/**
 * Write standard html5 index file
 * @param file_name output file name 
 */
 void writeIndexHtmlFile( const std::string & file_name )
 {
   std::ofstream file( file_name.c_str() ) ;
   if( ! file.good() )
   {
    std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
    exit( EXIT_FAILURE ) ; 
  }

  const char * const kHtmlHeader = 
  "<!DOCTYPE html>\n"
  "<html lang=\"en\">\n"
  "  <head>\n"
  "    <meta charset=\"utf-8\">\n"
  "    <title>WebGL viewer</title>\n"
  "    <link rel=\"stylesheet\" href=\"style/style.css\">\n"
  "\n"
  "    <!-- Point and line shaders -->\n"
  "    <script id=\"point_vShader\" type=\"x-shader/x-vertex\">\n"
  "    uniform mat4 uModelViewMatrix ;\n"
  "    uniform mat4 uProjectionMatrix ;\n"
  "    uniform mat4 uNormalMatrix;\n"
  "\n"
  "    attribute vec3 aNormal;\n"
  "    attribute vec3 aPosition;\n"
  "    attribute vec3 aColor;\n"
  "\n"
  "    varying vec4 vColor;\n"
  "\n" 
  "    void main( )\n"
  "    {\n"
  "      gl_PointSize = 1.0;\n"
  "      gl_Position = uProjectionMatrix * uModelViewMatrix * vec4( aPosition , 1.0 ) ;\n"
  "      vColor = vec4( aColor , 1.0 ) ;\n" 
  "    }\n"
  "    </script>\n"
  "\n"
  "    <script id=\"point_fShader\" type=\"x-shader/x-fragment\">\n"
  "    precision mediump float;\n"
  "    varying vec4 vColor ;\n"
  "\n"
  "    void main( )\n"
  "    {\n"
  "      gl_FragColor = vColor ;\n"
  "    }\n"
  "    </script>\n"
  "\n"
  "    <!-- Image (texture) shader -->\n"
  "    <script id=\"image_vShader\" type=\"x-shader/x-vertex\">\n"
  "    uniform mat4 uModelViewMatrix ;\n" 
  "    uniform mat4 uProjectionMatrix ;\n"
  "\n"
  "    attribute vec3 aPosition;\n"
  "    attribute vec2 aTexCoord;\n"
  "\n"
  "    varying highp vec2 vTexCoord;\n"
  "\n"
  "    void main( )\n"
  "    {\n"
  "      gl_Position = uProjectionMatrix * uModelViewMatrix * vec4( aPosition , 1.0 ) ;\n"
  "      vTexCoord = aTexCoord;\n"
  "    }\n"
  "    </script>\n"
  "\n"
  "    <script id=\"image_fShader\" type=\"x-shader/x-fragment\">\n"
  "    varying highp vec2 vTexCoord;\n" 
  "\n"
  "    uniform sampler2D uSampler ;\n" 
  "\n"
  "    void main( )\n"
  "    {\n"
  "      gl_FragColor = texture2D( uSampler , vec2(vTexCoord.s,vTexCoord.t));\n" 
  "    }\n" 
  "    </script>\n"
  "\n"
  "    <!-- other scripts -->\n"
  "    <script src=\"scripts/util.js\"></script>\n"
  "    <script src=\"scripts/matrix.js\"></script>\n"
  "    <script src=\"scripts/quaternion.js\"></script>\n"
  "    <script src=\"scripts/cloud.js\"></script>\n"
  "    <script src=\"scripts/shader.js\"></script>\n"
  "    <script src=\"scripts/cameras.js\"></script>\n"
  "    <script src=\"scripts/sceneData.js\"></script>\n"
  "    <script src=\"scripts/eventHandling.js\"></script>\n"
  "    <script src=\"scripts/main.js\"></script>\n"
  "  </head>\n" ;

  const char * const kHtmlBody = 
  "  <body onload=\"initGL()\">\n" 
  "    <header>\n"
  "      <h1>WebGL viewer</h1>\n" 
  "    </header>\n"
  "    <div id=\"container\">\n"
  "      <div id=\"divParams\">\n"
  "        <div class=\"divParamsHeader\">\n"
  "          <p>Display</p>\n"
  "        </div>\n"
  "        <div class=\"divParams\">\n"
  "          <div id=\"chkCams\" class=\"checkbox checked\">\n"
  "            <p>Cameras</p>\n"
  "          </div>\n"
  "          <div id=\"chkPcloud\" class=\"checkbox checked\">\n"
  "            <p>Point cloud</p>\n"
  "          </div>\n"
  "        </div>\n"
  "        <div class=\"divParamsFooter\"></div>\n"
  "        <div class=\"divParamsHeader\">\n"
  "          <p>Camera</p>\n"
  "        </div>\n"
  "        <div class=\"divParams\">\n"
  "          <div id=\"butResetCam\" class=\"button\">\n"
  "            <p>Reset</p>\n"
  "          </div>\n"
  "          <div id=\"chkAnimateCam\" class=\"checkbox checked\">\n"
  "            <p>Toggle animation</p>\n"
  "          </div>\n"
  "        </div>\n"
  "        <div class=\"divParamsFooter\"></div>\n"
  "      </div>\n"
  "\n"
  "      <canvas id=\"glCanvas\" width=\"1024\" height=\"768\">\n"
  "        Your browser does not support HTML5 canvas\n"
  "      </canvas>\n"
  "    </div>\n"
  "    <footer>\n"
  "      <p>OpenMVG WebGL viewer by Romuald PERROT</p>\n"
  "    </footer>\n"
  "  </body>\n"
  "</html>\n" ;

  file << kHtmlHeader ; 
  file << kHtmlBody ; 

  file.close() ; 
}

/**
* write event handling file 
* @param file_name output file name 
*/
void writeEventHandlingFile( const std::string & file_name )
{
  std::ofstream file( file_name.c_str() ) ;
  if( ! file.good() )
  {
    std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
    exit( EXIT_FAILURE ) ; 
  }

  const char * const kTogglePointCloudVisibleFunction =
  "function togglePointCloudVisible( )\n"
  "{\n"
  "  pCloud.ToggleVisible() ;\n"
  "\n"
  "  var elt = document.getElementById( \"chkPcloud\" ) ;\n"
  "  SwitchClass( elt , \"checked\" , \"unchecked\" ) ;\n"
  "}\n"
  "\n";

  const char * const kToggleCameraVisibleFunction =
  "function toggleCameraVisible( )\n"
  "{\n"
  "  for( var i = 0 ; i < cameras.length ; ++i )\n"
  "  {\n"
  "    cams[i].ToggleVisible() ; \n"
  "  }\n"
  "  var elt = document.getElementById( \"chkCams\" ) ;\n"
  "  SwitchClass( elt , \"checked\" , \"unchecked\" ) ;\n"
  "}\n"
  "\n";

  const char * const kToggleAnimation =
  "function ToggleAnimation()\n"
  "{\n"
  "  if( animateCam == true )\n"
  "  {\n"
  "    animateCam = false ;\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    animateCam = true ;\n"
  "  }\n"
  "  var elt = document.getElementById(\"chkAnimateCam\") ;\n"
  "  SwitchClass( elt , \"checked\" , \"unchecked\" ) ;\n"
  "}\n"
  "\n" ;

  const char * const kStopAnimation =
  "function StopAnimation()\n"
  "{\n"
  "  animateCam = false ; \n"
  "  var elt = document.getElementById(\"chkAnimateCam\") ;\n"
  "  if( HasClass( elt , \"checked\" ) )\n"
  "  {\n"
  "    RemoveClass( elt , \"checked\" ) ;\n"
  "    AddClass( elt , \"unchecked\" ) ;\n" 
  "  }\n"
  "}\n"
  "\n" ;

  const char * const kResetCamera =
  "function resetCamera()\n"
  "{\n"
  "  cameraModelView  = mat44.createLookAtMatrix( cam_center[0] , cam_center[1] , cam_center[2] , scene_center[0] , scene_center[1] , scene_center[2] , 0 , -1 , 0 );\n"
  "}\n"
  "\n" ;

  const char * const kKeyEvent =
  "function keyEvent( anEvent )\n"
  "{\n"
  "  var keyCode = anEvent.keyCode ;\n"
  "  var isValidKey = false ;\n"
  "\n"
  "  var camDir   = mat44.ExtractCameraForward( cameraModelView ) ; \n"
  "  var camRight = mat44.ExtractCameraRight( cameraModelView ) ;\n"
  "\n"
  "  if( keyCode == 90 ) // z\n"
  "  {\n"
  "    // Move forward \n"
  "    isValidKey = true ; \n"
  "\n"
  "    var tra = mat44.createTranslationMatrix( camDir[0] , camDir[1] , camDir[2] ) ; \n"
  "    cameraModelView = mat44.mul( tra , cameraModelView ) ; \n"
  "  }\n"
  "  if( keyCode == 81 ) // q\n"
  "  {\n"
  "    // move left\n"
  "    isValidKey = true ; \n"
  "\n"
  "    var tra = mat44.createTranslationMatrix( camRight[0] , camRight[1] , camRight[2] ) ; \n"
  "    cameraModelView = mat44.mul( tra , cameraModelView ) ;\n" 
  "  }\n"
  "  if( keyCode == 68 ) // d\n"
  "  {\n"
  "    // move right \n"
  "    isValidKey = true ; \n"
  "\n"
  "    var tra = mat44.createTranslationMatrix( -camRight[0] , -camRight[1] , -camRight[2] ) ; \n"
  "    cameraModelView = mat44.mul( tra , cameraModelView ) ; \n"
  "  }\n"
  "  if( keyCode == 83 ) // s\n"
  "  {\n"
  "    // move backward\n"
  "    isValidKey = true ; \n"
  "\n"
  "    var tra = mat44.createTranslationMatrix( -camDir[0] , -camDir[1] , -camDir[2] ) ; \n"
  "    cameraModelView = mat44.mul( tra , cameraModelView ) ; \n"
  "  }\n"
  "\n"
  "  if( isValidKey )\n"
  "  {\n"
  "    StopAnimation() ;\n"
  "  }\n"
  "}\n"
  "\n" ;

  const char * const kMouseVariables = 
  "/* Global variables handling mouse mouvement */\n"
  "var oldMouseX ;\n"
  "var oldMouseY ;\n"
  "var isMouseClicked = false ;\n"
  "\n" ;

  const char * const kMouseDown = 
  "/* Handle mouse click : start of dragging */\n"
  "function mouseDown( anEvent )\n"
  "{\n"
  "  var e = anEvent ? anEvent : window.event ; \n"
  "\n"
  "  oldMouseX = e.screenX ;\n"
  "  oldMouseY = e.screenY ;\n"
  "\n"
  "  isMouseClicked = true ; \n"
  "}\n"
  "\n";

  const char * const kMouseUp =
  "/* Handle mouse click : end of dragging */\n"
  "function mouseUp( anEvent )\n"
  "{\n"
  "  isMouseClicked = false ; \n"
  "}"
  "\n" ;

  const char * const kMouseOut =
  "/* Handle mouse moving outside gl canvas */\n"
  "function mouseOut( anEvent )\n"
  "{\n"
  "  isMouseClicked = false ; \n"
  "}\n"
  "\n" ;

  const char * const kMouseDrag = 
  "/* Main mouse gesture : change camera orientation */\n"
  "function mouseDrag( anEvent )\n"
  "{\n"
  "  if( isMouseClicked == false )\n"
  "    return ;\n"
  "\n"
  "  StopAnimation() ; \n"
  "\n"
  "  var e = anEvent ? anEvent : window.event ; \n"
  "\n"
  "  var deltaX = e.screenX - oldMouseX ;\n"
  "  var deltaY = e.screenY - oldMouseY ; \n"
  "\n"
  "  var camRight = mat44.ExtractCameraRight( cameraModelView ) ;\n"
  "  var camUp    = mat44.ExtractCameraUp( cameraModelView ) ; \n"
  "  var camPos   = mat44.ExtractCameraPos( cameraModelView ) ; \n"
  "\n"
  "  var tra     = mat44.createTranslationMatrix( -camPos[0] , - camPos[1] , - camPos[2] ) ;\n"
  "  var invtra  = mat44.createTranslationMatrix( camPos[0] , camPos[1] , camPos[2] ) ;\n"
  "\n"
  "  var q1 = quat.createRotation( camRight , -0.001 * deltaY ) ;\n"
  "  var q2 = quat.createRotation( camUp , -0.001 * deltaX ) ; \n"
  "  var q = quat.mul( q2 , q1 ) ; \n"
  "  var rot = quat.ToMatrix( q ) ; \n"
  "\n"
  "  var geom    = mat44.mul( tra , mat44.mul( rot , invtra ) ) ;\n"
  "\n"
  "  cameraModelView = mat44.mul( geom , cameraModelView ) ; \n"
  "\n"
  "  oldMouseX = e.screenX ;\n"
  "  oldMouseY = e.screenY ; \n"
  "}\n"
  "\n" ; 

  const char * const kSetupListeners =
  "function setupListeners()\n"
  "{\n"
  "  var chkCams   = document.getElementById( \"chkCams\" ) ;\n"
  "  var chkPcloud = document.getElementById( \"chkPcloud\" ) ;\n"
  "  var butResetCam = document.getElementById( \"butResetCam\" ) ;\n"
  "  var chkCamsAnim = document.getElementById( \"chkAnimateCam\" ) ;\n"
  "  var glWindow    = document.getElementById( \"glCanvas\" ) ; \n"
  "\n"
  "  window.addEventListener(\"keydown\", this.keyEvent , false ) ;\n"
  "\n"
  "  chkPcloud.onclick = togglePointCloudVisible ; \n"
  "  chkCams.onclick   = toggleCameraVisible ;\n"
  "  butResetCam.onclick = resetCamera ;\n"
  "  chkAnimateCam.onclick = ToggleAnimation ;\n"
  "\n"
  "  glWindow.onmousedown = mouseDown ;\n"
  "  glWindow.onmouseup = mouseUp ; \n"
  "  glWindow.onmousemove = mouseDrag ;\n"
  "  glWindow.onmouseout = mouseOut ; \n"
  "}\n"
  "\n";

  file << kTogglePointCloudVisibleFunction ;
  file << kToggleCameraVisibleFunction ;
  file << kResetCamera ; 
  file << kToggleAnimation ; 
  file << kStopAnimation ; 
  file << kKeyEvent ; 
  file << kMouseVariables ;
  file << kMouseDown ;
  file << kMouseUp ; 
  file << kMouseOut ;
  file << kMouseDown ;
  file << kMouseDrag ; 
  file << kSetupListeners ; 

  file.close() ; 
}

/**
 * write standard css file 
 * @param file_name output file name 
 */
 void writeCSSFile( const std::string & file_name )
 {
   std::ofstream file( file_name.c_str() ) ;
   if( ! file.good() )
   {
    std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
    exit( EXIT_FAILURE ) ; 
  }

  const char * const kCSSContent = 
  "*\n"
  "{\n"
  "  color: rgb(146,146,146) ;\n"
  "  margin: 0px ;\n"
  "  padding: 0px ;\n"
  "}\n"
  "\n"
  "body\n"
  "{\n"
  "  background-color: rgb(68,68,68) ;\n"
  "}\n"
  "\n"
  "#glCanvas\n"
  "{\n"
  "  display: block ;\n"
  "  width: 1024px ;\n"
  "  margin-left: auto ;\n"
  "  margin-right: auto ;\n"
  "  margin-top: 20px ;\n"
  "}\n" 
  "\n"
  "header\n"
  "{\n"
  "  background-color: rgb(40,40,40) ;\n"
  "  width: 100% ;\n"
  "  height: 60px ;\n"
  "  line-height: 60px ;\n"
  "  text-align: center ;\n"
  "  box-shadow: 0 0px 15px 15px rgb(40,40,40) ;\n"
  "}\n"
  "\n"
  "footer\n"
  "{\n"
  "  background-color: rgb(40,40,40) ;\n"
  "  position: fixed ;\n"
  "  bottom: 0px ;\n"
  "  width: 100% ;\n"
  "  height: 60px ;\n"
  "  line-height: 60px ;\n"
  "  text-align: center ;\n"
  "  box-shadow: 0 0px 15px 15px rgb(40,40,40) ;\n"
  "}\n" ;

  const char * const kCSSMenuParams =
  "#divParams\n"
  "{\n"
  "  position: absolute;\n"
  "  margin-left:10px;\n"
  "  margin-top:20px;\n"
  "  width:250px ;\n"
  "}\n"
  "\n"
  ".divParamsHeader\n"
  "{\n"
  "  background-color: rgb(40,40,40) ;\n"
  "  color: rgb(190,190,190);\n"
  "  text-align: center;\n"
  "  font-size: 1.5em;\n"
  "  border-top-left-radius: 15px;\n"
  "  border-top-right-radius: 15px;\n"
  "}\n"
  "\n"
  ".divParamsFooter\n"
  "{\n"
  "  background-color: rgb(40,40,40) ;\n"
  "  height: 0.7em;\n"
  "  margin-bottom: 10px;\n"
  "  border-bottom-right-radius: 15px;\n"
  "  border-bottom-left-radius: 15px;\n"
  "}\n"
  "\n"
  ".divParams\n"
  "{\n"
  "  border-left: 1px solid black;\n"
  "  border-right: 1px solid black ;\n"
  "  font-size: 1em;\n"
  "}\n"
  "\n"
  ".checkbox\n"
  "{\n"
  "  margin-right: 0.2em;\n"
  "  text-align: right;\n"
  "  cursor: pointer;\n"
  "  font-size: 1.1em;\n"
  "  background-color: rgb(68,68,68);\n"
  "\n"
  "  -webkit-touch-callout: none;\n"
  "  -webkit-user-select: none;\n"
  "  -khtml-user-select: none;\n"
  "  -moz-user-select: none;\n"
  "  -ms-user-select: none;\n"
  "  user-select: none;\n"
  "}\n"
  "\n"
  ".checkbox:hover\n"
  "{\n"
  "  background-color: rgb(80,80,80);\n"
  "}\n"
  "\n"
  ".checked p\n"
  "{\n"
  "  color:lightgreen;\n"
  "}\n"
  "\n"
  ".unchecked p\n"
  "{\n"
  "  color:red;\n"
  "}\n"
  "\n"
  ".button\n"
  "{\n"
  "  margin-right: 0.2em ; \n"
  "  text-align: right ;\n" 
  "  cursor: pointer;\n"
  "  font-size: 1.1em;\n"
  "  background-color: rgb(68,68,68);\n"
  "\n"
  "  -webkit-touch-callout: none;\n"
  "  -webkit-user-select: none;\n"
  "  -khtml-user-select: none;\n"
  "  -moz-user-select: none;\n"
  "  -ms-user-select: none;\n"
  "  user-select: none;\n"
  "}\n"
  "\n"
  ".button:hover p\n"
  "{\n"
  "  background-color: rgb(80,80,80);\n"
  "  font-weight: bold;\n"
  "  color:red;\n"
  "}\n"
  "\n";


  file << kCSSContent ;
  file << kCSSMenuParams ; 

  file.close() ;
}

/**
  * Write project scene specific data
  * @param doc the document containing scene data 
  * @param file_name output file name 
  */
void writeSceneDataFile( Document & doc , const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  // export sparse point cloud data 
	file << "cloud_data = new Float32Array( " << doc._vec_points.size() * 2 /* ie / 3 * 6 */ << ");" << std::endl;
	for( int i = 0 ; i < doc._vec_points.size() / 3 ; ++i )
	{
    // position 
		file << "cloud_data[" << 6 * i << "] = " << doc._vec_points[ 3 * i ] << ";" ;
		file << " cloud_data[" << 6 * i + 1 << "] = " << doc._vec_points[ 3 * i + 1 ] << ";" ;
		file << " cloud_data[" << 6 * i + 2 << "] = " << doc._vec_points[ 3 * i + 2 ] << ";" << std::endl ;
    // color (default : white)
		file << "cloud_data[" << 6 * i + 3 << "] = " << 1.0 << ";" ;
		file << " cloud_data[" << 6 * i + 4 << "] = " << 1.0 << ";" ;
		file << " cloud_data[" << 6 * i + 5 << "] = " << 1.0 << ";" << std::endl ;
	}

  // Compute camera distance in order to see point cloud 
	float center[3] = { 0 , 0 , 0 };
	float min_x =   std::numeric_limits<float>::max() ;
	float max_x = - std::numeric_limits<float>::max() ; 

	float min_y =   std::numeric_limits<float>::max() ;
	float max_y = - std::numeric_limits<float>::max() ; 

	for( int i = 0 ; i < doc._vec_points.size() / 3 ; ++i )
	{
		const float x = doc._vec_points[ 3 * i ] ; 
		min_x = std::min( min_x , x ) ; 
		max_x = std::max( max_x , x ) ; 

		const float y = doc._vec_points[ 3 * i + 1 ] ;
		min_y = std::min( min_y , y ) ; 
		max_y = std::max( max_y , y ) ;  

		center[0] += x ;
		center[1] += y ;
		center[2] += doc._vec_points[ 3 * i + 2 ] ; 
	}

	center[0] /= (float) doc._vec_points.size() ;
	center[1] /= (float) doc._vec_points.size() ; 
	center[2] /= (float) doc._vec_points.size() ;

	const float d_x = max_x - min_x ;
	const float d_y = max_y - min_y ;
	const float delta = std::max( d_x , d_y ) ;  
	const float dist  = delta / ( 2.0 * tan( D2R(60.0) / 2.0 ) ) ; // assume camera fov is 60 degrees

	const float center_cam[3] = { center[0] , center[1] , center[2] - dist } ; 

	file << "var scene_center = new Float32Array( 3 );" << std::endl ; 
	file << "scene_center[0] = " << center[0] << " ;" << std::endl ; 
	file << "scene_center[1] = " << center[1] << " ;" << std::endl ; 
	file << "scene_center[2] = " << center[2] << " ;" << std::endl ;	
	file << std::endl ; 
	file << "var cam_center  = new Float32Array( 3 );" << std::endl ;
	file << "cam_center[0] = " << center_cam[0] + 0.1 << " ;" << std::endl ; 
	file << "cam_center[1] = " << center_cam[1] + 0.1 << " ;" << std::endl ; 
	file << "cam_center[2] = " << center_cam[2] << " ;" << std::endl ;  

	file << "var cameraModelView  = mat44.createLookAtMatrix( cam_center[0] , cam_center[1] , cam_center[2] , scene_center[0] , scene_center[1] , scene_center[2] , 0 , -1 , 0 );" << std::endl ; 
	file << "var cameraProjection = mat44.createPerspectiveMatrix( 60.0 , 1024.0 / 768.0 , 0.1 , 1000.0 ) ;" << std::endl ; 

	// Export Cameras data
	file << "var cameras = new Array(" << doc._map_camera.size() << ") ;" << std::endl ; 
	for( int i = 0 ; i < doc._map_camera.size() ; ++i )
	{
		const Mat3 & R = doc._map_camera[i]._R ;
		const Mat3 & K = doc._map_camera[i]._K ; 
		const Vec3 & C = doc._map_camera[i]._C ; 

		file << "cameras[" << i << "] = new Object() ;" << std::endl ; 
		file << "cameras[" << i << "].position = new Float32Array( 3 ) ;" << std::endl ; 
    // position 
		file << "cameras[" << i << "].position[0] = " << C(0) << " ;" << std::endl ;
		file << "cameras[" << i << "].position[1] = " << C(1) << " ;" << std::endl ; 
		file << "cameras[" << i << "].position[2] = " << C(2) << " ;" << std::endl ; 
    // associated image 
		file << "cameras[" << i << "].imageName = \"./images/\" + (\"00000000\" +" << i << ").slice(-8) + \".jpg\"" << std::endl ; 
    // extremal points of iamge plane  
		file << "cameras[" << i << "].imagePlane = new Float32Array( 12 ) ;" << std::endl ; 
		const Vec3 p0 = R.transpose() * K.inverse() * Vec3( 0.0 , 0.0 , 1.0 ) + C ;
		const Vec3 p1 = R.transpose() * K.inverse() * Vec3( doc._map_imageSize[i].first , 0.0 , 1.0 ) + C ; 
		const Vec3 p2 = R.transpose() * K.inverse() * Vec3( doc._map_imageSize[i].first , doc._map_imageSize[i].second , 1.0 ) + C ; 
		const Vec3 p3 = R.transpose() * K.inverse() * Vec3( 0 , doc._map_imageSize[i].second , 1.0 ) + C ;
		file << "cameras[" << i << "].imagePlane[0] = " << p0(0) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[1] = " << p0(1) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[2] = " << p0(2) << " ; " << std::endl ;

		file << "cameras[" << i << "].imagePlane[3] = " << p1(0) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[4] = " << p1(1) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[5] = " << p1(2) << " ; " << std::endl ;
		
		file << "cameras[" << i << "].imagePlane[6] = " << p2(0) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[7] = " << p2(1) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[8] = " << p2(2) << " ; " << std::endl ;

		file << "cameras[" << i << "].imagePlane[9] = " << p3(0) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[10] = " << p3(1) << " ; " << std::endl ;
		file << "cameras[" << i << "].imagePlane[11] = " << p3(2) << " ; " << std::endl ;
	}


	file.close() ; 
}

/**
* Write minimal matrix class file 
* @param file_name output file name 
*/
void writeMatrixWebGLFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kStandardVectorOperations = 
  "/* L2 norm of a vector */\n"
  "function norm( x , y , z )\n"
  "{\n"
  "  return Math.sqrt( x * x + y * y + z * z ) ;\n" 
  "}\n"
  "\n"
  "/* dot product of two vectors */\n"
  "function dot( x1 , y1 , z1 , x2 , y2 , z2 )\n"
  "{\n;"
  "  return x1 * x2 + y1 * y2 + z1 * z2 ;"
  "}"
  "\n" ;
  
  const char * const kMatrixClass = 
  "/* Matrix class */\n"
  "var mat44 = {} ;\n"
  "\n";

  const char * const kLookAtMatrix = 
  "/* OpenGL like look at function */\n"
  "mat44.createLookAtMatrix = function( eyex , eyey , eyez , centerx , centery , centerz , _upx , _upy , _upz )\n"
  "{\n"
  "  var dstx = eyex - centerx ;\n"
  "  var dsty = eyey - centery ;\n"
  "  var dstz = eyez - centerz ;\n"
  "\n"
  "  var inv_norm = 1.0 / norm( dstx , dsty , dstz ) ;\n"
  "  dstx *= inv_norm ;\n"
  "  dsty *= inv_norm ;\n"
  "  dstz *= inv_norm ;\n"
  "\n"
  "  var upx = _upx ;\n"
  "  var upy = _upy ;\n"
  "  var upz = _upz ;\n"
  "\n"
  "  inv_norm = 1.0 / norm( upx , upy , upz ) ;\n"
  "  upx *= inv_norm ;\n"
  "  upy *= inv_norm ;\n"
  "  upz *= inv_norm ;\n"
  "\n"
  "  var rightx ;\n"
  "  var righty ;\n"
  "  var rightz ;\n"
  "\n"
  "  rightx = upy * dstz - dsty * upz ;\n"
  "  righty = upz * dstx - dstz * upx ;\n"
  "  rightz = upx * dsty - dstx * upy ;\n"
  "\n"
  "  inv_norm = 1.0 / norm( rightx , righty , rightz ) ;\n"
  "  rightx *= inv_norm ;\n"
  "  righty *= inv_norm ;\n"
  "  rightz *= inv_norm ;\n"
  "\n"
  "  upx = -righty * dstz + rightz * dsty ;\n"
  "  upy = -rightz * dstx + rightx * dstz ;\n"
  "  upz = -rightx * dsty + righty * dstx ;\n"
  "\n"
  "  var res = new Float32Array( 16 ) ;\n"
  "  res[0] = rightx ;\n"
  "  res[1] = upx ;\n"
  "  res[2] = dstx ;\n"
  "  res[3] = 0.0 ;\n"
  "  res[4] = righty ;\n"
  "  res[5] = upy ;\n"
  "  res[6] = dsty ;\n"
  "  res[7] = 0.0 ;\n"
  "  res[8] = rightz ;\n"
  "  res[9] = upz ;\n"
  "  res[10] = dstz ;\n"
  "  res[11] = 0.0 ;\n" 
  "  res[12] = - dot( rightx , righty , rightz , eyex , eyey , eyez ) ;\n"
  "  res[13] = - dot( upx , upy , upz , eyex , eyey , eyez ) ;\n"
  "  res[14] = - dot( dstx , dsty , dstz , eyex , eyey , eyez ) ;\n"
  "  res[15] = 1.0 ;\n" 
  "\n"
  "  return res ;\n"
  "};\n"
  "\n" ;

  const char * const kPerspectiveMatrix =
  "/* OpenGL (glu)perspective like matrix */\n"
  "mat44.createPerspectiveMatrix = function( fov , aspect , near , far )\n"
  "{\n"
  "  var range  = Math.tan( fov * Math.PI / 360.0 ) * near ;\n"
  "  var left   = -range * aspect ;\n"
  "  var right  = range * aspect ;\n"
  "  var bottom = - range ;\n"
  "  var top   = range ;\n"
  "\n"
  "  var out = new Float32Array( 16 ) ;\n"
  "\n"
  "  out[0] =  ( 2.0 * near ) / ( right - left ) ;\n"
  "  out[1] = 0.0 ;\n"
  "  out[2] = 0.0 ;\n"
  "  out[3] = 0.0 ;\n"
  "  out[4] = 0.0 ;\n"
  "  out[5] = ( 2.0 * near) / (top - bottom) ;\n"
  "  out[6] = 0.0 ;\n"
  "  out[7] = 0.0 ;\n"
  "  out[8] = (right + left) / ( right - left );\n"
  "  out[9] = (top + bottom) / ( top - bottom ) ;\n"
  "  out[10] = - (far + near ) / ( far - near ) ;\n"
  "  out[11] = -1.0 ;\n"
  "  out[12] = 0.0 ;\n"
  "  out[13] = 0.0 ;\n"
  "  out[14] = - ( 2.0 * far * near ) / ( far - near ) ;\n"
  "  out[15] = 0.0 ;\n"
  "\n"
  "  return out ;\n"
  "};\n"
  "\n" ;

  const char * const kMatrixMultiplication = 
  "/* Matrix multiplication */\n"
  "mat44.mul = function( a , b )\n"
  "{\n"
  "  var out = new Float32Array(16) ;\n"
  "  for( var i = 0 ; i < 4 ; ++i )\n" 
  "  {\n"
  "    for( var j = 0 ; j < 4 ; ++j )\n"
  "    {\n"
  "      var idx = i * 4 + j ;\n"
  "      out[idx] = 0.0 ;\n"
  "\n"
  "        for( var k = 0 ; k < 4 ; ++k )\n"
  "        {\n"
  "          out[idx] += a[i*4+k] * b[k*4+j];\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "\n"
  "  return out ;\n"
  "};\n"
  "\n" ;

  const char * const kMatrixInversion = 
  "/* Matrix inversion */\n"
  "mat44.invert = function( a )\n"
  "{\n"
  "  /* see intel paper: streaming SIMD Extensions - Inverse of 4x4 matrix   */\n"
  "  var  tmp = new Float32Array( 16 ) ; /* temp array for pairs             */\n"
  "  var  src = new Float32Array( 16 ) ; /* array of transpose source matrix */\n"
  "  var out = new Float32Array( 16 ) ;\n"
  "\n"
  "  /* transpose matrix */\n"
  "  for (var i = 0; i < 4; ++i)\n"
  "  {\n"
  "    src[i]        = a[i*4];\n"
  "    src[i + 4]    = a[i*4 + 1];\n"
  "    src[i + 8]    = a[i*4 + 2];\n"
  "    src[i + 12]   = a[i*4 + 3];\n"
  "  }\n"
  "\n"
  "/* calculate pairs for first 8 elements (cofactors) */\n"
  "  tmp[0]  = src[10] * src[15];\n"
  "  tmp[1]  = src[11] * src[14];\n"
  "  tmp[2]  = src[9]  * src[15];\n"
  "  tmp[3]  = src[11] * src[13];\n"
  "  tmp[4]  = src[9]  * src[14];\n"
  "  tmp[5]  = src[10] * src[13];\n"
  "  tmp[6]  = src[8]  * src[15];\n"
  "  tmp[7]  = src[11] * src[12];\n"
  "  tmp[8]  = src[8]  * src[14];\n"
  "  tmp[9]  = src[10] * src[12];\n"
  "  tmp[10] = src[8]  * src[13];\n"
  "  tmp[11] = src[9]  * src[12];\n"
  "  /* calculate first 8 elements (cofactors) */\n"
  "  out[0]  = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];\n"
  "  out[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];\n"
  "  out[1]  = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];\n"
  "  out[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];\n"
  "  out[2]  = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];\n"
  "  out[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];\n"
  "  out[3]  = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];\n"
  "  out[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];\n"
  "  out[4]  = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];\n"
  "  out[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];\n"
  "  out[5]  = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];\n"
  "  out[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];\n"
  "  out[6]  = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];\n"
  "  out[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];\n"
  "  out[7]  = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];\n"
  "  out[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];\n"
  "  /* calculate pairs for second 8 elements (cofactors) */\n"
  "  tmp[0]  = src[2]*src[7];\n"
  "  tmp[1]  = src[3]*src[6];\n"
  "  tmp[2]  = src[1]*src[7];\n"
  "  tmp[3]  = src[3]*src[5];\n"
  "  tmp[4]  = src[1]*src[6];\n"
  "  tmp[5]  = src[2]*src[5];\n"
  "  tmp[6]  = src[0]*src[7];\n"
  "  tmp[7]  = src[3]*src[4];\n"
  "  tmp[8]  = src[0]*src[6];\n"
  "  tmp[9]  = src[2]*src[4];\n"
  "  tmp[10] = src[0]*src[5];\n"
  "  tmp[11] = src[1]*src[4];\n"
  "  /* calculate second 8 elements (cofactors) */\n"
  "  out[8]  = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];\n"
  "  out[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];\n"
  "  out[9]  = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];\n"
  "  out[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];\n"
  "  out[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];\n"
  "  out[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];\n"
  "  out[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];\n"
  "  out[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];\n"
  "  out[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];\n"
  "  out[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];\n"
  "  out[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];\n"
  "  out[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];\n"
  "  out[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];\n"
  "  out[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];\n"
  "  out[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];\n"
  "  out[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];\n"
  "  /* calculate determinant */\n"
  "  var det=src[0]*out[0]+src[1]*out[1]+src[2]*out[2]+src[3]*out[3];\n"
  "  /* calculate matrix inverse */\n"
  "  det = 1 / det;\n"
  "  for (var j = 0; j < 16; ++j )\n"
  "   out[j] *= det;\n"
  "\n"
  "  return out ;\n"
  "};\n"
  "\n" ;

  const char * const kTranslationMatrix =
  "/* Translation matrix */\n"
  "mat44.createTranslationMatrix = function( dx , dy , dz )\n"
  "{\n"
  "  var out = new Float32Array( 16 ) ;\n"
  "\n"
  "  out[0] = 1.0 ;\n"
  "  out[1] = 0.0 ;\n"
  "  out[2] = 0.0 ;\n"
  "  out[3] = 0.0 ;\n"
  "\n"
  "  out[4] = 0.0 ;\n"
  "  out[5] = 1.0 ;\n"
  "  out[6] = 0.0 ;\n"
  "  out[7] = 0.0 ;\n"
  "\n"
  "  out[8] = 0.0 ;\n"
  "  out[9] = 0.0 ;\n"
  "  out[10] = 1.0 ;\n"
  "  out[11] = 0.0 ;\n"
  "\n"
  "  out[12] = dx ;\n"
  "  out[13] = dy ;\n"
  "  out[14] = dz ;\n"
  "  out[15] = 1.0 ;\n"
  "\n"
  "  return out ;\n"
  "};\n"
  "\n" ;

  const char * const kRotationYMatrix = 
  "/* Rotation about y axis - angle_rad is rotation angle in radian*/\n"
  "mat44.createRotateYMatrix = function( angle_rad )\n"
  "{\n"
  "  var s = Math.sin( angle_rad ) ;\n"
  "  var c = Math.cos( angle_rad ) ;\n"
  "  var out = new Float32Array( 16 ) ;\n"
  "\n"
  "  out[0] = c ;\n"
  "  out[1] = 0 ;\n"
  "  out[2] = s ;\n"
  "  out[3] = 0 ;\n"
  "\n"
  "  out[4] = 0.0 ;\n"
  "  out[5] = 1.0 ;\n"
  "  out[6] = 0.0 ;\n"
  "  out[7] = 0.0 ;\n"
  "\n"
  "  out[8] = -s ;\n"
  "  out[9] = 0.0 ;\n"
  "  out[10] = c ;\n"
  "  out[11] = 0.0 ;\n"
  "\n"
  "  out[12] = 0.0 ;\n"
  "  out[13] = 0.0 ;\n"
  "  out[14] = 0.0 ;\n"
  "  out[15] = 1.0 ;\n"
  "\n"
  "  return out ;\n"
  "};\n"
  "\n" ;


  const char * const kExtractCameraPos = 
  "mat44.ExtractCameraPos = function( aMatrix )\n"
  "{\n"
  "  var inve = mat44.invert( aMatrix ) ;\n"
  "  var res = new Float32Array( 3 ) ;\n"
  "  res[0] = inve[12] ;\n"
  "  res[1] = inve[13] ;\n"
  "  res[2] = inve[14] ;\n" 
  "  return res ;\n" 
  "}\n"
  "\n" ;

  const char * const kExtractCameraUp = 
  "mat44.ExtractCameraUp = function( aMatrix )\n"
  "{\n"
  "  var res = new Float32Array( 3 ) ;\n"
  "  res[0] = aMatrix[1] ;\n"
  "  res[1] = aMatrix[5] ;\n"
  "  res[2] = aMatrix[9] ;\n"
  "  return res ;\n" 
  "}\n"
  "\n";

  const char * const kExtractCameraForward = 
  "mat44.ExtractCameraForward = function( aMatrix )\n"
  "{\n"
  "  var res = new Float32Array( 3 ) ;\n"
  "  res[0] = aMatrix[2] ;\n"
  "  res[1] = aMatrix[6] ;\n"
  "  res[2] = aMatrix[10] ;\n"
  "  return res ;\n"
  "}\n"
  "\n";

  const char * const kExtractCameraRight = 
  "mat44.ExtractCameraRight = function( aMatrix )\n"
  "{\n"
  "  var res = new Float32Array( 3 ) ;\n"
  "  res[0] = aMatrix[0] ;\n"
  "  res[1] = aMatrix[4] ;\n"
  "  res[2] = aMatrix[8] ;\n"
  "  return res ;\n"
  "}\n"
  "\n" ;

  file << kStandardVectorOperations ;
  file << kMatrixClass ;
  file << kLookAtMatrix ; 
  file << kPerspectiveMatrix ; 
  file << kMatrixMultiplication ; 
  file << kMatrixInversion ; 
  file << kTranslationMatrix ; 
  file << kRotationYMatrix ; 
  file << kExtractCameraPos ;
  file << kExtractCameraUp ;
  file << kExtractCameraForward ;
  file << kExtractCameraRight ; 

  file.close() ; 
}

/**
 * write quaternion class 
 * @param file_name output file name  
 */
void writeQuaternionFile( const std::string & file_name )
{
  std::ofstream file( file_name.c_str() ) ;
  if( ! file.good() )
  {
    std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
    exit( EXIT_FAILURE ) ; 
  }

  const char * const kQuaternionConstructor = 
  "/* minimum quaternion class */\n"
  "var quat = {} ;\n"
  "\n" ;

  const char * const kQuaternionCreateRotation =
  "/* Create a rotation quaternion from an axis and an angle */\n"
  "quat.createRotation = function( anAxis , anAngle )\n"
  "{\n"
  "  var inv_norm = 1.0 / norm( anAxis[0] , anAxis[1] , anAxis[2] ) ;\n"
  "  var axis = new Float32Array( 3 ) ;\n"
  "  axis[0] = anAxis[0] * inv_norm ;\n"
  "  axis[1] = anAxis[1] * inv_norm ;\n"
  "  axis[2] = anAxis[2] * inv_norm ;\n"
  "\n"
  "  var theta = anAngle / 2.0 ;\n"
  "  var c = Math.cos( theta ) ;\n"
  "  var s = Math.sin( theta ) ;\n"
  "\n"
  "  var res = new Float32Array( 4 ) ;\n"
  "  res[0] = c ;\n"
  "  res[1] = s * axis[0] ;\n"
  "  res[2] = s * axis[1] ;\n"
  "  res[3] = s * axis[2] ;\n"
  "  return res ; \n"
  "}\n"
  "\n" ;

  const char * const kQuaternionMultiplication =
  "/* Quaternion rotation */\n"
  "quat.mul = function( aQuat1 , aQuat2 )\n"
  "{\n"
  "  var a = aQuat1[0] ;\n"
  "  var b = aQuat1[1] ;\n"
  "  var c = aQuat1[2] ;\n"
  "  var d = aQuat1[3] ;\n"
  "\n"
  "  var ap = aQuat2[0] ;\n"
  "  var bp = aQuat2[1] ;\n"
  "  var cp = aQuat2[2] ;\n"
  "  var dp = aQuat2[3] ;\n"
  "\n"
  "  var res = new Float32Array( 4 ) ;\n"
  "  res[0] = a * ap - b * bp - c * cp - d * dp ;\n"
  "  res[1] = a * bp + b * ap + c * dp - d * cp ;\n"
  "  res[2] = a * cp - b * dp + c * ap + d * bp ;\n"
  "  res[3] = a * dp + b * cp - c * bp + d * ap ; \n"
  "\n"
  "  return res ; \n"
  "}\n"
  "\n" ;

  const char * const kQuaternionToMatrix =
  "/* Convert quaternion to rotation matrix */\n"
  "quat.ToMatrix = function( aQuaternion )\n"
  "{\n"
  "  var a = aQuaternion[0] ;\n"
  "  var b = aQuaternion[1] ;\n"
  "  var c = aQuaternion[2] ;\n"
  "  var d = aQuaternion[3] ;\n"
  "\n"
  "  var a2 = a * a ;\n"
  "  var b2 = b * b ;\n"
  "  var c2 = c * c ;\n"
  "  var d2 = d * d ; \n"
  "\n"
  "  var ab = a * b ;\n"
  "  var ac = a * c ; \n"
  "  var ad = a * d ; \n"
  "  var bc = b * c ; \n"
  "  var bd = b * d ;\n"
  "  var cd = c * d ; \n"
  "\n"
  "  var res = new Float32Array( 16 ) ;\n"
  "\n"
  "  res[0] = a2 + b2 - c2 - d2 ;\n"
  "  res[1] = 2.0 * ( bc - ad ) ;\n"
  "  res[2] = 2.0 * ( bd + ac ) ;\n"
  "  res[3] = 0.0 ;\n"
  "\n"
  "  res[4] = 2.0 * ( bc + ad ) ;\n"
  "  res[5] = a2 - b2 + c2 - d2 ;\n"
  "  res[6] = 2.0 * ( cd - ab ) ;\n"
  "  res[7] = 0.0 ;\n"
  "\n"
  "  res[8] = 2.0 * ( bd - ac ) ;\n"
  "  res[9] = 2.0 * ( cd + ab ) ;\n"
  "  res[10] = a2 - b2 - c2 + d2 ;\n"
  "  res[11] = 0.0 ;\n"
  "\n"
  "  res[12] = 0.0 ;\n"
  "  res[13] = 0.0 ;\n"
  "  res[14] = 0.0 ;\n"
  "  res[15] = 1.0 ;\n"
  "\n"
  "  return res ;\n"   
  "}\n"
  "\n" ;

  file << kQuaternionConstructor ;
  file << kQuaternionMultiplication ;
  file << kQuaternionCreateRotation ;
  file << kQuaternionToMatrix ; 

  file.close() ; 
}



/** 
  * Write main WebGL script code 
  * @param file_name output file name 
  */
void writeMainWebGLFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kGlobalVariables = 
  "/* WebGL Context */\n"
  "var gl ;\n"
  "\n"
  "/* Shaders */\n"
  "var pointShaderProgram ;\n"
  "var surfaceShaderProgram ;\n"
  "\n"
  "/* Scene objects : point cloud and cameras */\n"
  "var pCloud ;\n"
  "var cams ;\n"
  "\n"
  "/* Animation */\n"
  "var animateCam ;\n"
  "\n";

  const char * const kWebGLInitialization = 
  "/* WebGL init */\n"
  "function initGL()\n"
  "{\n" 
  "  var canvas = document.getElementById(\"glCanvas\") ;\n"
  "  gl         = canvas.getContext(\"webgl\") || canvas.getContext(\"experimental-webgl\") ;\n"
  "  if( ! gl )\n"
  "  {\n"
  "    alert(\"could not load WebGL context\");\n"
  "    return ;\n"
  "  }\n"
  "  setup( ) ; /* setup scene */\n"
  "  setupListeners() ; /* setup interface event */\n"
  "  update( ) ; /* launch main loop */\n"
  "}\n"
  "\n" ;

  const char * const kSetupWebGL = 
  "/* setup scene */\n"
  "function setup( )\n"
  "{\n"
  "  gl.clearColor( 0.0 , 0.0 , 0.0 , 1.0 ) ;\n"
  "  gl.enable( gl.DEPTH_TEST ) ;\n"
//  "  gl.enable( gl.VERTEX_PROGRAM_POINT_SIZE);\n"
  "  gl.depthFunc( gl.LEQUAL ) ;\n"
  "\n"
  "  pointShaderProgram = new Shader( gl , \"point_vShader\" , \"point_fShader\" ) ;\n"
  "  pCloud = new PointCloud( gl , pointShaderProgram , cloud_data ) ;\n"
  "  surfaceShaderProgram = new Shader( gl , \"image_vShader\" , \"image_fShader\" ) ;\n"
  "  cams   = new Array( cameras.length ) ;\n"
  "  for( var i = 0 ; i < cameras.length ; ++i )\n"
  "  {\n"
  "    cams[i] = new Camera( gl , surfaceShaderProgram , pointShaderProgram , cameras[i].position , cameras[i].imagePlane , cameras[i].imageName ) ;\n"
  "  }\n"
  "  animateCam = true ; \n"
  "}\n"
  "\n" ;

  const char * const kMainLoop = 
  "/* Main loop : animate then render scene */\n"
  "function update()\n"
  "{\n"
  "  requestAnimFrame( update );\n"
  "\n"
  "  /* Rotate about center of scene */\n"
  "  if( animateCam )\n"
  "  {\n"
  "    var tra     = mat44.createTranslationMatrix( - scene_center[0] , - scene_center[1] , - scene_center[2] ) ;\n"
  "    var inv_tra = mat44.createTranslationMatrix( scene_center[0] , scene_center[1] , scene_center[2] );\n"
  "    var rot     = mat44.createRotateYMatrix( 0.005 ) ;\n"
  "    var geom    = mat44.mul( tra , mat44.mul( rot , inv_tra ) );\n"
  "    cameraModelView = mat44.mul( geom , cameraModelView ) ;\n"
  "  }\n"
  "\n"
  "  render() ;\n"
  "}\n"
  "\n" ;

  const char * const kRenderFunction = 
  "/* Rendering function */\n"
  "function render( )\n"
  "{\n" 
  "  gl.clear( gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT ) ;\n"
  "\n"
  "  pointShaderProgram.Enable( gl ) ;\n"
  "  var mv_loc = gl.getUniformLocation( pointShaderProgram.shad , \"uModelViewMatrix\" ) ;\n"
  "  var pj_loc = gl.getUniformLocation( pointShaderProgram.shad , \"uProjectionMatrix\" ) ;\n"
  "  gl.uniformMatrix4fv( mv_loc , false , cameraModelView ) ;\n"
  "  gl.uniformMatrix4fv( pj_loc , false , cameraProjection ) ;\n"
  "  pCloud.render( gl ) ;\n"
  "  for( var i = 0 ; i < cameras.length ; ++i )\n"
  "  {\n"
  "    cams[i].renderLines( gl );\n"
  "  }\n"
  "  surfaceShaderProgram.Enable( gl ) ; \n"
  "  var mv_loc = gl.getUniformLocation( surfaceShaderProgram.shad , \"uModelViewMatrix\" ) ;\n"
  "  var pj_loc = gl.getUniformLocation( surfaceShaderProgram.shad , \"uProjectionMatrix\" ) ;\n"
  "  gl.uniformMatrix4fv( mv_loc , false , cameraModelView ) ;\n"
  "  gl.uniformMatrix4fv( pj_loc , false , cameraProjection ) ;\n"
  "  for( var i = 0 ; i < cameras.length ; ++i )\n"
  "  {\n"
  "    cams[i].render( gl ) ;\n"
  "  }\n"
  "}\n"
  "\n" ;

  file << kGlobalVariables ; 
  file << kWebGLInitialization ; 
  file << kSetupWebGL ; 
  file << kMainLoop ; 
  file << kRenderFunction ; 

	file.close() ;
}


/** 
 * Write camera object class file 
 * @param file_name output file name 
 */
void WriteCameraFile( const std::string & file_name ) 
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kCameraConstructor = 
  "function Camera( aGLContext , aShader , aShaderLine , aPosition , aPlanePos , aFileName )\n"
  "{\n" 
  "  this.pos           = aPosition ;\n"
  "  this.planePos      = aPlanePos ;\n"
  "  this.shader        = aShader ;\n"
  "  this.shaderLines   = aShaderLine;\n"
  "  this.imageFileName = aFileName ;\n"
  "\n"
  "  /* Build Image plane geometry */\n"
  "  // array = 2 triangles : 0 1 2 - 0 2 3\n"
  "  var buff_img = new Float32Array( 30 ) ;\n"
  "\n"
  "  // 0\n"
  "  buff_img[0] = aPlanePos[0] ;\n"
  "  buff_img[1] = aPlanePos[1] ;\n"
  "  buff_img[2] = aPlanePos[2] ;\n"
  "  buff_img[3] = 0.0 ;\n"
  "  buff_img[4] = 0.0 ;\n"
  "\n"
  "  // 1\n"
  "  buff_img[5] = aPlanePos[3] ;\n"
  "  buff_img[6] = aPlanePos[4] ;\n"
  "  buff_img[7] = aPlanePos[5] ;\n"
  "  buff_img[8] = 1.0 ;\n"
  "  buff_img[9] = 0.0 ;\n"
  "\n"
  "  // 2\n"
  "  buff_img[10] = aPlanePos[6] ;\n"
  "  buff_img[11] = aPlanePos[7] ;\n"
  "  buff_img[12] = aPlanePos[8] ;\n"
  "  buff_img[13] = 1.0 ;\n"
  "  buff_img[14] = 1.0 ;\n"
  "\n"
  "  // 0\n"
  "  buff_img[15] = aPlanePos[0] ;\n"
  "  buff_img[16] = aPlanePos[1] ;\n"
  "  buff_img[17] = aPlanePos[2] ;\n"
  "  buff_img[18] = 0.0 ;\n"
  "  buff_img[19] = 0.0 ;\n"
  "\n"
  "  // 2\n"
  "  buff_img[20] = aPlanePos[6] ;\n" 
  "  buff_img[21] = aPlanePos[7] ;\n"
  "  buff_img[22] = aPlanePos[8] ;\n"
  "  buff_img[23] = 1.0 ;\n"
  "  buff_img[24] = 1.0 ;\n"
  "\n"
  "  // 3\n"
  "  buff_img[25] = aPlanePos[9] ;\n"
  "  buff_img[26] = aPlanePos[10] ;\n"
  "  buff_img[27] = aPlanePos[11] ;\n"
  "  buff_img[28] = 0.0 ;\n"
  "  buff_img[29] = 1.0 ;\n"
  "\n"
  "  var buff_lines = new Float32Array( 48 ) ;\n"
  "  buff_lines[0] = aPosition[0] ;\n"
  "  buff_lines[1] = aPosition[1] ;\n"
  "  buff_lines[2] = aPosition[2] ;\n"
  "  buff_lines[3] = 1.0 ;\n"
  "  buff_lines[4] = 1.0 ;\n"
  "  buff_lines[5] = 1.0 ;\n"
  "  buff_lines[6] = aPlanePos[0] ;\n"
  "  buff_lines[7] = aPlanePos[1] ;\n"
  "  buff_lines[8] = aPlanePos[2] ;\n"
  "  buff_lines[9] = 1.0 ;\n"
  "  buff_lines[10] = 1.0 ;\n"
  "  buff_lines[11] = 1.0 ;\n"
  "  buff_lines[12] = aPosition[0] ;\n"
  "  buff_lines[13] = aPosition[1] ;\n"
  "  buff_lines[14] = aPosition[2] ;\n"
  "  buff_lines[15] = 1.0 ;\n"
  "  buff_lines[16] = 1.0 ;\n"
  "  buff_lines[17] = 1.0 ;\n"
  "  buff_lines[18] = aPlanePos[3] ;\n"
  "  buff_lines[19] = aPlanePos[4] ;\n"
  "  buff_lines[20] = aPlanePos[5] ;\n"
  "  buff_lines[21] = 1.0 ;\n"
  "  buff_lines[22] = 1.0 ;\n"
  "  buff_lines[23] = 1.0 ;\n"
  "  buff_lines[24] = aPosition[0] ;\n"
  "  buff_lines[25] = aPosition[1] ;\n"
  "  buff_lines[26] = aPosition[2] ;\n"
  "  buff_lines[27] = 1.0 ;\n"
  "  buff_lines[28] = 1.0 ;\n"
  "  buff_lines[29] = 1.0 ;\n"
  "  buff_lines[30] = aPlanePos[6] ;\n"
  "  buff_lines[31] = aPlanePos[7] ;\n"
  "  buff_lines[32] = aPlanePos[8] ;\n"
  "  buff_lines[33] = 1.0 ;\n"
  "  buff_lines[34] = 1.0 ;\n"
  "  buff_lines[35] = 1.0 ;\n"
  "  buff_lines[36] = aPosition[0] ;\n"
  "  buff_lines[37] = aPosition[1] ;\n"
  "  buff_lines[38] = aPosition[2] ;\n"
  "  buff_lines[39] = 1.0 ;\n"
  "  buff_lines[40] = 1.0 ;\n"
  "  buff_lines[41] = 1.0 ;\n"
  "  buff_lines[42] = aPlanePos[9] ;\n"
  "  buff_lines[43] = aPlanePos[10] ;\n"
  "  buff_lines[44] = aPlanePos[11] ;\n"
  "  buff_lines[45] = 1.0 ;\n"
  "  buff_lines[46] = 1.0 ;\n"
  "  buff_lines[47] = 1.0 ;\n"
  "  // create vbo\n"
  "  this.vbo = aGLContext.createBuffer() ;\n"
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;\n"
  "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , buff_img , aGLContext.STATIC_DRAW ) ;\n"
  "\n"
  "  this.vboLines = aGLContext.createBuffer() ;\n"
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vboLines ) ;\n"
  "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , buff_lines , aGLContext.STATIC_DRAW ) ;\n"
  "  // create texture file\n"
  "  var tex       = gl.createTexture() ;\n"
  "  var img = new Image() ;\n"
  "  this.tex = tex ;\n"
  "  this.img = img ;\n"
  "  img.onload    = function() { handleTexture( aGLContext , img , tex ); } ;"
  "  img.src       = aFileName ;\n"
  "  this.nbelt     = 6 ;\n"
  "  this.visible   = true ;\n"
  "}\n"
  "\n";

  const char * const kToggleVisible =
  "Camera.prototype.ToggleVisible = function()\n"
  "{\n"
  "  if( this.visible == false )\n"
  "    this.visible = true ;\n"
  "  else\n"
  "    this.visible = false ;\n"
  "}\n"
  "\n";

  const char * const kTextureHandling =
  "/* texture handling */\n"
  "function handleTexture( aGLContext , anImage , aTexture )\n" 
  "{\n" 
  "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, aTexture);\n"
  "  aGLContext.texImage2D(aGLContext.TEXTURE_2D, 0, aGLContext.RGBA, aGLContext.RGBA, aGLContext.UNSIGNED_BYTE, anImage);\n"
  "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_MAG_FILTER, aGLContext.LINEAR);\n"
  "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_MIN_FILTER, aGLContext.LINEAR);\n"
  "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_WRAP_S, aGLContext.CLAMP_TO_EDGE);\n"
  "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_WRAP_T, aGLContext.CLAMP_TO_EDGE);\n"
  "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, null);\n"
  "}\n"
  "\n" ;

  const char * const kRenderImage = 
  "/* Main rendering function (ie render image plane) */\n"
  "Camera.prototype.render = function( aGLContext )\n" 
  "{\n"
  "  if( this.visible == false )\n"
  "    return ;\n"
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;\n"
  "\n"
  "  aGLContext.enableVertexAttribArray( this.shader.attribPos ) ;\n"
  "  aGLContext.enableVertexAttribArray( this.shader.attribTex ) ;\n"
  "\n"
  "  aGLContext.vertexAttribPointer( this.shader.attribPos , 3 , aGLContext.FLOAT , false , 20 , 0 ) ;\n"
  "  aGLContext.vertexAttribPointer( this.shader.attribTex , 2 , aGLContext.FLOAT , false , 20 , 12 ) ;\n"
  "\n"
  "  aGLContext.activeTexture(aGLContext.TEXTURE0);\n"
  "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, this.tex);\n"
  "  aGLContext.uniform1i(aGLContext.getUniformLocation(this.shader.shad, \"uSampler\"), 0);\n"
  "\n"
  "  aGLContext.drawArrays( aGLContext.TRIANGLES , 0 , this.nbelt ) ;\n"
  "}\n"
  "\n";

  const char * const kRenderLines = 
  "/* Render camera structure (ie support lines)*/\n"
  "Camera.prototype.renderLines = function( aGLContext )\n"
  "{\n"
  "  if( this.visible == false )\n"
  "    return ;\n"
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vboLines ) ;\n"
  "\n"
  "  aGLContext.enableVertexAttribArray( this.shaderLines.attribPos ) ;\n"
  "  aGLContext.enableVertexAttribArray( this.shaderLines.attribCol ) ;\n"
  "\n"
  "  aGLContext.vertexAttribPointer( this.shaderLines.attribPos , 3 , aGLContext.FLOAT , false , 24 , 0 ) ;\n"
  "  aGLContext.vertexAttribPointer( this.shaderLines.attribCol , 3 , aGLContext.FLOAT , false , 24 , 12 ) ;\n"
  "\n"
  "  aGLContext.drawArrays( aGLContext.LINES , 0 , 8 ) ;\n"
  "}\n"
  "\n" ;

  file << kCameraConstructor ;
  file << kToggleVisible ; 
  file << kTextureHandling ; 
  file << kRenderImage ; 
  file << kRenderLines ; 

  file.close() ; 
}

/**
* Write shader file 
* @param file_name output file name 
*/
void writeShaderFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kShaderConstructor = 
  "function Shader( aGLContext , aVertexID , aFragID )\n"  
  "{\n" 
  "  this.shad = setupShader( aGLContext , aVertexID , aFragID ) ;\n" 
  "\n"
  "  this.attribPos = aGLContext.getAttribLocation( this.shad , \"aPosition\" ) ;\n"
  "  this.attribTex = aGLContext.getAttribLocation( this.shad , \"aTexCoord\" ) ;\n"
  "  this.attribNor = aGLContext.getAttribLocation( this.shad , \"aNormal\" ) ;\n"
  "  this.attribCol = aGLContext.getAttribLocation( this.shad , \"aColor\" ) ;\n"
  "}\n"
  "\n";

  const char * const kShaderEnable = 
  "/* Activate shader */\n"
  "Shader.prototype.Enable = function( aGLContext )\n" 
  "{\n"
  "  aGLContext.useProgram( this.shad ) ;\n"
  "}\n"
  "\n" ;

  const char * const kShaderDisable = 
  "/* deactivate shader */\n"
  "Shader.prototype.Disable = function( aGLContext )\n"
  "{\n"
  "  aGLContext.useProgram( 0 ) ;\n"
  "}\n"
  "\n" ;

  const char * const kCompileShader =
  "/* shader compilation (vertex of fragment)*/\n"
  "function compileShader( aGLContext , aShaderSource , aShaderType )\n"
  "{\n"  
  "  var shad = aGLContext.createShader( aShaderType ) ;\n"
  "  aGLContext.shaderSource( shad , aShaderSource ) ;\n"
  "  aGLContext.compileShader( shad ) ;\n"
  "\n"
  "  var ok = aGLContext.getShaderParameter( shad , aGLContext.COMPILE_STATUS ) ;\n"
  "  if( ! ok )\n"
  "  {\n"
  "    throw \"could not compile shader : \" + aGLContext.getShaderInfoLog( shad ) ;\n"
  "  }\n"
  "\n"
  "  return shad ;\n"
  "}\n"
  "\n";

  const char * const kSetupShader = 
  "/* Setup shader program (compiles both vertex and fragment and link) */\n"
  "function setupShader( aGLContext , aVertexID , aFragID )\n"  
  "{\n" 
  "  var vElt = document.getElementById( aVertexID ) ;\n"
  "  if( ! vElt )\n"
  "  {\n" 
  "    throw \"could not find vertex shader\" + aVertexID ;\n" 
  "  }\n" 
  "  var vStr  = vElt.text ;\n" 
  "  var vShad = compileShader( aGLContext , vStr , aGLContext.VERTEX_SHADER ) ;\n" 
  "  var fElt = document.getElementById( aFragID ) ;\n"
  "  if( ! fElt )\n"
  "  {\n" 
  "    throw \"could not find fragment shader \" + aFragID ;\n"
  "  }\n" 
  "  var fStr  = fElt.text ;\n" 
  "  var fShad = compileShader( aGLContext , fStr , aGLContext.FRAGMENT_SHADER ) ;\n"
  "  var pgm = aGLContext.createProgram( ) ;\n"
  "  aGLContext.attachShader( pgm , vShad ) ;\n"
  "  aGLContext.attachShader( pgm , fShad ) ;\n"
  "  aGLContext.linkProgram( pgm ) ;\n" 
  "  var ok = aGLContext.getProgramParameter( pgm , aGLContext.LINK_STATUS ) ;\n"
  "  if( ! ok )\n"
  "  {\n" 
  "    throw \"pgm failed to link\" + aGLContext.getProgramInfoLog( pgm ) ;\n"
  "  }\n"
  "  return pgm ;\n"
  "}\n" 
  "\n" ;

  file << kShaderConstructor ; 
  file << kShaderEnable ; 
  file << kShaderDisable ; 
  file << kCompileShader ; 
  file << kSetupShader ; 

  file.close() ; 
}

/**
  * Write Point cloud class file 
  * @param file_name output file name 
  */
void writePointCloudFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kPointCloudConstructor = 
  "function PointCloud( aGLContext , aShader , aPointData )\n" 
  "{\n" 
  "  this.pointData = aPointData ;\n"
  "  this.shader    = aShader ;\n" 
  "  this.nbelt     = aPointData.length / 6 ; /* 6 = 3 position + 3 color*/\n"
  "  this.visible   = true ;\n"
  "\n"
  "  this.vbo = aGLContext.createBuffer( ) ;\n" 
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;\n"
  "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , aPointData , aGLContext.STATIC_DRAW ) ;\n"
  "}\n"
  "\n" ;

  const char * const kToggleVisible = 
  "PointCloud.prototype.ToggleVisible = function( )\n"
  "{\n"
  "  if( this.visible == false )\n"
  "    this.visible = true ;\n" 
  "  else\n"
  "    this.visible = false ;\n"
  "}\n"
  "\n";

  const char * const kRenderFunction = 
  "/* Render point cloud */\n"
  "PointCloud.prototype.render = function( aGLContext )\n"  
  "{\n"
  "  if( this.visible == false )\n"
  "    return ;\n"
  "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;\n" 
  "\n"
  "  aGLContext.enableVertexAttribArray( this.shader.attribPos ) ;\n"
  "  aGLContext.enableVertexAttribArray( this.shader.attribCol ) ;\n" 
  "\n"
  "  aGLContext.vertexAttribPointer( this.shader.attribPos , 3 , aGLContext.FLOAT , false , 24 , 0 ) ;\n"
  "  aGLContext.vertexAttribPointer( this.shader.attribCol , 3 , aGLContext.FLOAT , false , 24 , 12 ) ;\n" 
  "\n"
  "  aGLContext.drawArrays( aGLContext.POINTS , 0 , this.nbelt ) ;\n" 
  "}\n"
  "\n" ; 

  file << kPointCloudConstructor ; 
  file << kToggleVisible ; 
  file << kRenderFunction ; 

	file.close() ; 
}


/** 
 * Write javascript util file 
 * @param file_name output file name 
 */
void writeUtilFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

  const char * const kRequestAnimFunction = 
  "/* Portable requestAnimFrame fonction */\n"
  "requestAnimFrame = (function()\n" 
  "{\n" 
  "  return window.requestAnimationFrame ||\n"
  "   window.webkitRequestAnimationFrame ||\n"
  "   window.mozRequestAnimationFrame ||\n"
  "   window.oRequestAnimationFrame ||\n"
  "   window.msRequestAnimationFrame ||\n" 
  "   function(/* function FrameRequestCallback */ callback, /* DOMElement Element */ element) {\n"
  "   window.setTimeout(callback, 1000/60);\n"
  " };\n"
  "})();\n"
  "\n" ;

  const char * const kHasClassFunction =
  "/* Check if DOM element has a class property */\n"
  "function HasClass( aNode , aClassName )\n"
  "{\n"
  "  if( aNode.className )\n"
  "  {\n"
  "    return aNode.className.match( new RegExp( '(\\\\s|^)' + aClassName + '(\\\\s|$)') ) ;\n"
  "  }\n"
  "  else\n"
  "  {\n"
  "    return false ;\n"
  "  }\n"
  "}\n"
  "\n" ;

  const char * const kAddClassFunction =
  "/* Add a classname to a DOM element */\n"
  "function AddClass( aNode , aClassName )\n"
  "{\n"
  "  aNode.className += \" \" + aClassName ;\n"
  "}\n"
  "\n" ;

  const char * const kRemoveClassFunction = 
  "/* Remove a classname to a DOM element */\n"
  "function RemoveClass( aNode , aClassName )\n"
  "{\n"
  "  if( HasClass( aNode , aClassName ) )\n"
  "  {\n"
  "    var reg = new RegExp('(\\\\s|^)' + aClassName + '(\\\\s|$)');\n"
  "    aNode.className = aNode.className.replace(reg, ' ');\n"
  "  }\n"
  "}\n"
  "\n" ;

  const char * const kSwitchClassFunction = 
  "/* Switch two class in a DOM element */\n"
  "function SwitchClass( aNode , aValue1 , aValue2 )\n"
  "{\n"
  "  if( HasClass( aNode , aValue1 ) )\n"
  "  {\n"
  "    RemoveClass( aNode , aValue1 ) ;\n"
  "    AddClass( aNode , aValue2 ) ;\n"
  "  }\n"
  "  else if( HasClass( aNode , aValue2 ) )\n"
  "  {\n"
  "    RemoveClass( aNode , aValue2 ) ;\n"
  "    AddClass( aNode , aValue1 ) ;\n"
  "  }\n"
  "}\n"
  "\n" ;

  file << kRequestAnimFunction ;  
  file << kHasClassFunction ;
  file << kAddClassFunction ;
  file << kRemoveClassFunction ;
  file << kSwitchClassFunction ;
  file.close() ; 
}

/**
 * Export images to image folder 
 * @param doc openSFM document 
 * @param sSfMDir openSFM output folder (for locating input images)
 * @param folder_image output location 
 */
void writeImagesFiles( Document & doc , const std::string & sSfMDir , const std::string & folder_image )
{
	// Initial Image folder 
	const std::string sImagePath = stlplus::folder_append_separator(sSfMDir) + "images" ;

	Image<RGBColor> image ;
	Image<RGBColor> ressamp ;
	for( int i = 0 ; i < doc._vec_imageNames.size() ; ++i )
	{
		const std::string & sImageName = doc._vec_imageNames[i];
   ReadImage( stlplus::create_filespec( sImagePath, sImageName).c_str(), &image );

	    // Need to redimension image to max 2048x2048 

   std::ostringstream os;
   os << std::setw(8) << std::setfill('0') << i ;


   int nw = 1024 ;
   int nh = 1024 ;
   RessamplingSize( nw , nh , doc._map_imageSize[i].first , doc._map_imageSize[i].second ) ;
   Ressample( image , nw , nh , ressamp ) ; 

   const std::string & file_name = stlplus::folder_append_separator( folder_image ) + os.str() + ".jpg" ; 
   WriteImage( file_name.c_str() , ressamp ) ;
 }
}


/**
 * Export project to a fully functionnal WebGL site
 * @param doc openSFM document containing project to export 
 * @param sSfMDir openSFM project directory 
 * @param outFolder output directory 
 */
void exportProjectToWebGL( Document & doc , const std::string & sSfMDir , const std::string & outFolder )
{
	// Create output's directory structure 
	prepareFolders( outFolder ) ; 

	// write html and css file 
	const std::string html_file_name = stlplus::folder_append_separator( outFolder ) + "index.html" ; 
	writeIndexHtmlFile( html_file_name ) ; 
	const std::string css_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "style" ) + "style.css" ; 
	writeCSSFile( css_file_name ) ;

	// write util js file 
	const std::string util_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "util.js" ; 
	writeUtilFile( util_js_file_name ) ; 

	// write matrix js file 
	const std::string matrix_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "matrix.js" ; 
	writeMatrixWebGLFile( matrix_js_file_name ) ; 

  const std::string quaternion_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "quaternion.js" ; 
  writeQuaternionFile( quaternion_js_file_name ) ; 

	// write shader js file
	const std::string shader_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "shader.js" ; 
	writeShaderFile( shader_js_file_name ) ; 

	// write camera webgl file
	const std::string camera_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "cameras.js" ; 
	WriteCameraFile( camera_js_file_name ) ; 

	// write pointcloud webgl file 
	const std::string pcloud_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "cloud.js" ; 
	writePointCloudFile( pcloud_js_file_name ) ; 

  // write event handling 
  const std::string event_handling_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "eventHandling.js" ; 
  writeEventHandlingFile( event_handling_file_name );

  // write scene data
	const std::string scene_data_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "sceneData.js" ; 
	writeSceneDataFile( doc , scene_data_file_name ) ; 

	// write main webgl file 
	const std::string main_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "main.js" ; 
	writeMainWebGLFile( main_js_file_name ) ; 

	// write images files
	const std::string image_folder = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "images" );
	writeImagesFiles( doc , sSfMDir , image_folder ) ;

}

int main( int argc, char ** argv )
{
	CmdLine cmd ;
	// Sfm directory (ie: output of openMVG)
	std::string sSfmDir ;
	// output of densification process 
	std::string sOutDir ; 

	cmd.add( make_option('i' , sSfmDir , "sfmdir") ) ;
	cmd.add( make_option('o' , sOutDir , "output") ) ;


  try 
  {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  }
  catch(const std::string& s)
  {
    std::cerr << "Usage : " << argv[0] << " " 
    << "[-i|--sfmdir, the SfM_output path] " 
    << "[-o|--output, the output directory] " 
    << std::endl ; 

    std::cerr << s << std::endl ; 

    return EXIT_FAILURE ; 
  }

  Document m_doc ;

  if( m_doc.load( sSfmDir ) )
  {
    exportProjectToWebGL( m_doc , sSfmDir , sOutDir ) ; 
  }
  else
  {
    return EXIT_FAILURE ; 
  }

  return EXIT_SUCCESS ; 
}