#include <iostream>
#include <cstdlib>
#include <iomanip>


#include "software/SfMViewer/document.h"
#include "third_party/cmdLine/cmdLine.h"
#include "openMVG/image/image.hpp"


/**
  *
  * Prepare WebGL folder 
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


// given a new_width and a new_height 
// recompute new_width and new_height in order to conserve initial ratio.  
static inline void RessamplingSize( int & new_width , int & new_height , const int oldWidth , const int oldHeight )
{
	if( oldWidth > oldHeight )
	{
		const double ratio = (double) new_width / (double) oldWidth ; 
		new_height = ratio * (double) oldHeight ; 
	}
	else
	{
		const double ratio = (double) new_height / (double) oldHeight ;
		new_width = ratio * (double) oldWidth ; 
	}
}

/**
 * Write standard html5 index file  
 */
void writeIndexHtmlFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	// header 
	file << "<!DOCTYPE html>" << std::endl ; 
	file << "<html lang=\"en\">" << std::endl ; 
	file << "  <head>" << std::endl ; 
	file << "    <meta charset=\"utf-8\">" << std::endl ; 
	file << "    <title>WebGL viewer</title>" << std::endl ; 
	file << "    <link rel=\"stylesheet\" href=\"style/style.css\">" << std::endl ; 

	// shaders 
	file << "    <script id=\"point_vShader\" type=\"x-shader/x-vertex\">" << std::endl ; 
	file << "    uniform mat4 uModelViewMatrix ;" << std::endl ; 
	file << "    uniform mat4 uProjectionMatrix ;" << std::endl ; 
	file << "    uniform mat4 uNormalMatrix;" << std::endl ; 
	file << std::endl ; 
	file << "    attribute vec3 aNormal;" << std::endl ; 
	file << "    attribute vec3 aPosition;" << std::endl ; 
	file << "    attribute vec3 aColor;" << std::endl ; 
	file << std::endl ; 
	file << "    varying vec4 vColor;" << std::endl ; 
	file << std::endl ; 
	file << "    void main( )" << std::endl;
	file << "    {" << std::endl;
	file << "      gl_Position = uProjectionMatrix * uModelViewMatrix * vec4( aPosition , 1.0 ) ;" << std::endl ;
	file << "      vColor = vec4( aColor , 1.0 ) ;" << std::endl ; 
	file << "    }" << std::endl ; 
	file << "    </script>" << std::endl ; 

	file << "    <script id=\"point_fShader\" type=\"x-shader/x-fragment\">" << std::endl ; 
	file << "    precision mediump float;" << std::endl ; 
	file << "    varying vec4 vColor ; " << std::endl ; 
	file << std::endl ; 
	file << "    void main( )" << std::endl ; 
	file << "    {" << std::endl ; 
	file << "      gl_FragColor = vColor ;" << std::endl ; 
	file << "    }" << std::endl ; 
	file << "    </script>" << std::endl ; 

	file << "    <script id=\"image_vShader\" type=\"x-shader/x-vertex\">" << std::endl ; 
	file << "    uniform mat4 uModelViewMatrix ;" << std::endl ; 
	file << "    uniform mat4 uProjectionMatrix ;" << std::endl ; 
	file << std::endl ; 
	file << "    attribute vec3 aPosition;" << std::endl ; 
	file << "    attribute vec2 aTexCoord;" << std::endl ; 
	file << std::endl ; 
	file << "    varying highp vec2 vTexCoord;" << std::endl ; 
	file << std::endl ; 
	file << "    void main( )" << std::endl ;
	file << "    {" << std::endl ; 
	file << "      gl_Position = uProjectionMatrix * uModelViewMatrix * vec4( aPosition , 1.0 ) ;" << std::endl ; 
	file << "      vTexCoord = aTexCoord;" << std::endl;
	file << "    }" << std::endl ; 
	file << "    </script>" << std::endl ; 

	file << "    <script id=\"image_fShader\" type=\"x-shader/x-fragment\">" << std::endl ; 
	file << "    varying highp vec2 vTexCoord;" << std::endl ; 
	file << std::endl ; 
	file << "    uniform sampler2D uSampler ;" << std::endl ; 
	file << std::endl ; 
	file << "    void main( )" << std::endl ; 
	file << "    {" << std::endl ; 
	file << "      gl_FragColor = texture2D( uSampler , vec2(vTexCoord.s,vTexCoord.t));" << std::endl ; 
	file << "    }" << std::endl ; 
	file << "    </script>" << std::endl ; 

	file << "    <script src=\"scripts/util.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/matrix.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/cloud.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/shader.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/cameras.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/sceneData.js\"></script>" << std::endl ; 
	file << "    <script src=\"scripts/main.js\"></script>" << std::endl ; 
	file << "  </head>" << std::endl ; 

	// body 
	file << "  <body onload=\"initGL()\">" << std::endl ; 
	file << "    <header>" << std::endl ; 
	file << "      <h1>WebGL viewer</h1>" << std::endl ; 
	file << "    </header>" << std::endl ; 
	file << std::endl ; 
	file << "    <canvas id=\"glCanvas\" width=\"1024\" height=\"768\">" << std::endl ; 
	file << "      Your browser does not support HTML5 canvas" << std::endl ; 
	file << "    </canvas>" << std::endl ; 
	file << std::endl ; 
	file << "    <footer>" << std::endl ; 
	file << "      <p>OpenMVG WebGL viewer by Romuald PERROT</p>" << std::endl ;
	file << "    </footer>" << std::endl ; 
	file << "  </body>" << std::endl ; 

	file << "</html>" << std::endl ; 

	file.close() ; 
}

/**
 * write standard css file 
 */
void writeCSSFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "*" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  color: rgb(146,146,146) ;" << std::endl ;
	file << "  margin: 0px ;" << std::endl ;
	file << "  padding: 0px ;" << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ;
	file << "body" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  background-color: rgb(68,68,68) ;" << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "#glCanvas" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  display: block ;" << std::endl ; 
	file << "  width: 1024px ;" << std::endl ; 
	file << "  margin-left: auto ;" << std::endl ; 
	file << "  margin-right: auto ;" << std::endl ; 
	file << "  margin-top: 20px ;" << std::endl ;  
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "header" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  background-color: rgb(40,40,40) ;" << std::endl ; 
	file << "  width: 100% ;" << std::endl ;
	file << "  height: 60px ;" << std::endl ; 
	file << "  line-height: 60px ;" << std::endl ; 
	file << "  text-align: center ;" << std::endl ; 
	file << "  box-shadow: 0 0px 15px 15px rgb(40,40,40) ;" << std::endl ;
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "footer" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  background-color: rgb(40,40,40) ;" << std::endl ; 
	file << "  position: fixed ;" << std::endl ; 
	file << "  bottom: 0px ;" << std::endl ; 
	file << "  width: 100% ;" << std::endl ; 
	file << "  height: 60px ;" << std::endl ; 
	file << "  line-height: 60px ;" << std::endl ; 
	file << "  text-align: center ;" << std::endl ; 
	file << "  box-shadow: 0 0px 15px 15px rgb(40,40,40) ;" << std::endl ; 
	file << "}" << std::endl ; 

	file.close() ; 
}

void writeSceneDataFile( Document & doc , const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "cloud_data = new Float32Array( " << doc._vec_points.size() * 2 /* ie / 3 * 6 */ << ");" << std::endl;
	for( int i = 0 ; i < doc._vec_points.size() / 3 ; ++i )
	{
		file << "cloud_data[" << 6 * i << "] = " << doc._vec_points[ 3 * i ] << ";" ;
		file << " cloud_data[" << 6 * i + 1 << "] = " << doc._vec_points[ 3 * i + 1 ] << ";" ;
		file << " cloud_data[" << 6 * i + 2 << "] = " << doc._vec_points[ 3 * i + 2 ] << ";" << std::endl ;
		file << "cloud_data[" << 6 * i + 3 << "] = " << 1.0 << ";" ;
		file << " cloud_data[" << 6 * i + 4 << "] = " << 1.0 << ";" ;
		file << " cloud_data[" << 6 * i + 5 << "] = " << 1.0 << ";" << std::endl ;
	}

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
	const float dist  = delta / ( 2.0 * tan( D2R(60.0) / 2.0 ) ) ;

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

	file << "var cameraModelView  = mat44.createLookAtMatrix( cam_center[0] , cam_center[1] , cam_center[2] , scene_center[0] , scene_center[1] , scene_center[2] , 0 , 1 , 0 );" << std::endl ; 
	file << "var cameraProjection = mat44.createPerspectiveMatrix( 40.0 , 1024.0 / 768.0 , 0.1 , 1000.0 ) ;" << std::endl ; 

	// Export Cameras
	file << "var cameras = new Array(" << doc._map_camera.size() << ") ;" << std::endl ; 
	for( int i = 0 ; i < doc._map_camera.size() ; ++i )
	{
		const Mat3 & R = doc._map_camera[i]._R ;
		const Mat3 & K = doc._map_camera[i]._K ; 
		const Vec3 & C = doc._map_camera[i]._C ; 

		file << "cameras[" << i << "] = new Object() ;" << std::endl ; 
		file << "cameras[" << i << "].position = new Float32Array( 3 ) ;" << std::endl ; 
		file << "cameras[" << i << "].position[0] = " << C(0) << " ;" << std::endl ;
		file << "cameras[" << i << "].position[1] = " << C(1) << " ;" << std::endl ; 
		file << "cameras[" << i << "].position[2] = " << C(2) << " ;" << std::endl ; 
		file << "cameras[" << i << "].imageName = \"./images/\" + (\"00000000\" +" << i << ").slice(-8) + \".jpg\"" << std::endl ; 
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

void writeMatrixWebGLFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "function norm( x , y , z )" << std::endl ;
	file << "{" << std::endl ;
	file << "  return Math.sqrt( x * x + y * y + z * z ) ;" << std::endl;
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "function dot( x1 , y1 , z1 , x2 , y2 , z2 )" << std::endl ;
	file << "{" << std::endl ;
	file << "  return x1 * x2 + y1 * y2 + z1 * z2 ;" << std::endl ; 
	file << "}" << std::endl ;
	file << std::endl ; 
	file << "var mat44 = {} ; " << std::endl ; 
	file << std::endl ; 

	// look at 
	file << "mat44.createLookAtMatrix = function( eyex , eyey , eyez , centerx , centery , centerz , _upx , _upy , _upz )" << std::endl ;
	file << "{" << std::endl ; 
	file << "  var dstx = eyex - centerx ;" << std::endl ;
	file << "  var dsty = eyey - centery ;" << std::endl ; 
	file << "  var dstz = eyez - centerz ;" << std::endl ; 
	file << std::endl ; 
	file << "  var inv_norm = 1.0 / norm( dstx , dsty , dstz ) ;" << std::endl;
	file << "  dstx *= inv_norm ;" << std::endl ; 
	file << "  dsty *= inv_norm ;" << std::endl ; 
	file << "  dstz *= inv_norm ;" << std::endl ; 
	file << std::endl ;
	file << "  var upx = _upx ;" << std::endl ; 
	file << "  var upy = _upy ;" << std::endl ; 
	file << "  var upz = _upz ;" << std::endl ; 
	file << std::endl ; 
	file << "  inv_norm = 1.0 / norm( upx , upy , upz ) ;" << std::endl ; 
	file << "  upx *= inv_norm ;" << std::endl ; 
	file << "  upy *= inv_norm ;" << std::endl ; 
	file << "  upz *= inv_norm ;" << std::endl ;
	file << std::endl ; 
	file << "  var rightx ;" << std::endl ; 
	file << "  var righty ;" << std::endl ; 
	file << "  var rightz ;" << std::endl ; 
	file << std::endl ; 
	file << "  rightx = upy * dstz - dsty * upz ;" << std::endl ; 
	file << "  righty = upz * dstx - dstz * upx ;" << std::endl ; 
	file << "  rightz = upx * dsty - dstx * upy ;" << std::endl ;
	file << std::endl ; 
	file << "  inv_norm = 1.0 / norm( rightx , righty , rightz ) ;" << std::endl ; 
	file << "  rightx *= inv_norm ;" << std::endl ;
	file << "  righty *= inv_norm ;" << std::endl ; 
	file << "  rightz *= inv_norm ;" << std::endl ; 
	file << std::endl ; 
	file << "  upx = righty * dstz - rightz * dsty ;" << std::endl ;
	file << "  upy = rightz * dstx - rightx * dstz ;" << std::endl ; 
	file << "  upz = rightx * dsty - righty * dstx ;" << std::endl ; 
	file << std::endl ; 
	file << "  var res = new Float32Array( 16 ) ;" << std::endl ; 
	file << "  res[0] = rightx ;" << std::endl ; 
	file << "  res[1] = upx ;" << std::endl ; 
	file << "  res[2] = dstx ;" << std::endl ; 
	file << "  res[3] = 0.0 ;" << std::endl ; 
	file << "  res[4] = righty ;" << std::endl ; 
	file << "  res[5] = upy ;" << std::endl ; 
	file << "  res[6] = dsty ;" << std::endl ; 
	file << "  res[7] = 0.0 ;" << std::endl ; 
	file << "  res[8] = rightz ;" << std::endl ; 
	file << "  res[9] = upz ;" << std::endl ; 
	file << "  res[10] = dstz ;" << std::endl ;
	file << "  res[11] = 0.0 ;" << std::endl ; 
	file << "  res[12] = - dot( rightx , righty , rightz , eyex , eyey , eyez ) ;" << std::endl ; 
	file << "  res[13] = - dot( upx , upy , upz , eyex , eyey , eyez ) ;" << std::endl ; 
	file << "  res[14] = - dot( dstx , dsty , dstz , eyex , eyey , eyez ) ;" << std::endl ;
	file << "  res[15] = 1.0 ;" << std::endl; 
	file << std::endl ;
	file << "  return res ;" << std::endl;
	file << "};" << std::endl ;
	file << std::endl ; 

	// perspective 
	file << "mat44.createPerspectiveMatrix = function( fov , aspect , near , far )" << std::endl ; 
	file << "{" << std::endl ;
	file << "  var range  = Math.tan( fov * Math.PI / 360.0 ) * near ;" << std::endl ; 
	file << "  var left   = -range * aspect ;" << std::endl ; 
	file << "  var right  = range * aspect ;" << std::endl ; 
	file << "  var bottom = - range ;" << std::endl ; 
	file << "  var top   = range ;" << std::endl ; 
	file << std::endl ; 
	file << "  var out = new Float32Array( 16 ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  out[0] =  ( 2.0 * near ) / ( right - left ) ;" << std::endl ; 
	file << "  out[1] = 0.0 ;" << std::endl ; 
	file << "  out[2] = 0.0 ;" << std::endl ; 
	file << "  out[3] = 0.0 ;" << std::endl ; 
	file << "  out[4] = 0.0 ;" << std::endl ; 
	file << "  out[5] = ( 2.0 * near) / (top- bottom) ;" << std::endl ; 
	file << "  out[6] = 0.0 ;" << std::endl ; 
	file << "  out[7] = 0.0 ;" << std::endl ; 
	file << "  out[8] = (right + left) / ( right - left );" << std::endl ; 
	file << "  out[9] = (top + bottom) / ( top - bottom ) ;" << std::endl ; 
	file << "  out[10] = - (far + near ) / ( far - near ) ;" << std::endl ; 
	file << "  out[11] = -1.0 ;" << std::endl ; 
	file << "  out[12] = 0.0 ;" << std::endl ; 
	file << "  out[13] = 0.0 ;" << std::endl ; 
	file << "  out[14] = - ( 2.0 * far * near ) / ( far - near ) ;" << std::endl ; 
	file << "  out[15] = 0.0 ;" << std::endl ; 
	file << std::endl ; 
	file << "  return out ;" << std::endl ; 
	file << "};" << std::endl ; 

	// multiplication 
	file << "mat44.mul = function( a , b )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  var out = new Float32Array(16) ;" << std::endl ; 
	file << "  for( var i = 0 ; i < 4 ; ++i )" << std::endl ; 
	file << "  {" << std::endl ; 
	file << "    for( var j = 0 ; j < 4 ; ++j )" << std::endl ; 
	file << "    {" << std::endl ; 
	file << "      var idx = i * 4 + j ;" << std::endl ; 
	file << "      out[idx] = 0.0 ;" << std::endl ; 
	file << std::endl ; 
	file << "        for( var k = 0 ; k < 4 ; ++k )" << std::endl ; 
	file << "        {" << std::endl ; 
	file << "          out[idx] += a[i*4+k] * b[k*4+j];" << std::endl ; 
	file << "        }" << std::endl ; 
	file << "      }" << std::endl ; 
	file << "    }" << std::endl ; 
	file << std::endl ; 
	file << "  return out ;" << std::endl ; 
	file << "};" << std::endl ;

	// inverse (really still needed  ?)
	file << "mat44.invert = function( a )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  // /* intel : streaming SIMD Extensions - Inverse of 4x4 matrix */" << std::endl ; 
	file << "  var  tmp = new Float32Array( 16 ) ; /* temp array for pairs                      */ " << std::endl ; 
	file << "  var  src = new Float32Array( 16 ) ; /* array of transpose source matrix */" << std::endl ; 
	file << "  var out = new Float32Array( 16 ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  /* transpose matrix */" << std::endl ; 
	file << "  for (var i = 0; i < 4; ++i) " << std::endl ;
	file << "  {" << std::endl ; 
	file << "    src[i]        = a[i*4]; " << std::endl ; 
	file << "    src[i + 4]    = a[i*4 + 1];" << std::endl ; 
	file << "    src[i + 8]    = a[i*4 + 2];" << std::endl ; 
	file << "    src[i + 12]   = a[i*4 + 3];" << std::endl ; 
	file << "  }" << std::endl ; 
	file << std::endl ; 
	file << "/* calculate pairs for first 8 elements (cofactors) */" << std::endl ; 
	file << "  tmp[0]  = src[10] * src[15];" << std::endl ; 
	file << "  tmp[1]  = src[11] * src[14];" << std::endl ; 
	file << "  tmp[2]  = src[9]  * src[15];" << std::endl ; 
	file << "  tmp[3]  = src[11] * src[13];" << std::endl ;
	file << "  tmp[4]  = src[9]  * src[14];" << std::endl ; 
	file << "  tmp[5]  = src[10] * src[13];" << std::endl ; 
	file << "  tmp[6]  = src[8]  * src[15];" << std::endl ;  
	file << "  tmp[7]  = src[11] * src[12];" << std::endl ;
	file << "  tmp[8]  = src[8]  * src[14];" << std::endl ; 
	file << "  tmp[9]  = src[10] * src[12];" << std::endl ; 
	file << "  tmp[10] = src[8]  * src[13];" << std::endl ; 
	file << "  tmp[11] = src[9]  * src[12];" << std::endl ; 
	file << "  /* calculate first 8 elements (cofactors) */" << std::endl ; 
	file << "  out[0]  = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];" << std::endl ; 
	file << "  out[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];" << std::endl ; 
	file << "  out[1]  = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];" << std::endl ; 
	file << "  out[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];" << std::endl ;
	file << "  out[2]  = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];" << std::endl ; 
	file << "  out[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];" << std::endl ; 
	file << "  out[3]  = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];" << std::endl ; 
	file << "  out[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];" << std::endl ; 
	file << "  out[4]  = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];" << std::endl ; 
	file << "  out[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];" << std::endl ; 
	file << "  out[5]  = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];" << std::endl ; 
	file << "  out[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];" << std::endl ; 
	file << "  out[6]  = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];" << std::endl ; 
	file << "  out[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];" << std::endl ; 
	file << "  out[7]  = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];" << std::endl ; 
	file << "  out[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];" << std::endl ; 
	file << "  /* calculate pairs for second 8 elements (cofactors) */" << std::endl ;
	file << "  tmp[0]  = src[2]*src[7];" << std::endl ; 
	file << "  tmp[1]  = src[3]*src[6];" << std::endl ;
	file << "  tmp[2]  = src[1]*src[7];" << std::endl ; 
	file << "  tmp[3]  = src[3]*src[5];" << std::endl ; 
	file << "  tmp[4]  = src[1]*src[6];" << std::endl ; 
	file << "  tmp[5]  = src[2]*src[5];" << std::endl ;
	file << "  tmp[6]  = src[0]*src[7];" << std::endl ; 
	file << "  tmp[7]  = src[3]*src[4];" << std::endl ; 
	file << "  tmp[8]  = src[0]*src[6];" << std::endl ; 
	file << "  tmp[9]  = src[2]*src[4];" << std::endl ; 
	file << "  tmp[10] = src[0]*src[5];" << std::endl ; 
	file << "  tmp[11] = src[1]*src[4];" << std::endl ;
	file << "  /* calculate second 8 elements (cofactors) */" << std::endl ; 
	file << "  out[8]  = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; " << std::endl ;
	file << "  out[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15]; " << std::endl ;
	file << "  out[9]  = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; " << std::endl ;
	file << "  out[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; " << std::endl ;
	file << "  out[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15]; " << std::endl ;
	file << "  out[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; " << std::endl ;
	file << "  out[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14]; " << std::endl ;
	file << "  out[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; " << std::endl ;
	file << "  out[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9]; " << std::endl ;
	file << "  out[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; " << std::endl ;
	file << "  out[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; " << std::endl ;
	file << "  out[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; " << std::endl ;
	file << "  out[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; " << std::endl ;
	file << "  out[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; " << std::endl ;
	file << "  out[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; " << std::endl ;
	file << "  out[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8]; " << std::endl ;
	file << "  /* calculate determinant */ " << std::endl ;
	file << "  var det=src[0]*out[0]+src[1]*out[1]+src[2]*out[2]+src[3]*out[3]; " << std::endl ;
	file << "  /* calculate matrix inverse */ " << std::endl ;
	file << "  det = 1 / det; " << std::endl ;
	file << "  for (var j = 0; j < 16; ++j ) " << std::endl ;
	file << "  	out[j] *= det; " << std::endl ;
	file << std::endl; 
	file << "  return out ; " << std::endl ;
    file << "};" << std::endl ;

    // translation
    file << "mat44.createTranslationMatrix = function( dx , dy , dz )" << std::endl ; 
    file << "{" << std::endl ; 
    file << "  var out = new Float32Array( 16 ) ; " << std::endl ; 
    file << std::endl ; 
    file << "  out[0] = 1 ;" << std::endl ; 
    file << "  out[1] = 0 ;" << std::endl ; 
    file << "  out[2] = 0 ; " << std::endl ; 
    file << "  out[3] = 0 ; " << std::endl ; 
    file << std::endl ; 
    file << "  out[4] = 0 ;" << std::endl ; 
    file << "  out[5] = 1 ;" << std::endl ; 
    file << "  out[6] = 0 ;" << std::endl ; 
    file << "  out[7] = 0 ;" << std::endl ; 
    file << std::endl ; 
    file << "  out[8] = 0 ;" << std::endl ; 
    file << "  out[9] = 0 ;" << std::endl ; 
    file << "  out[10] = 1 ;" << std::endl ; 
    file << "  out[11] = 0 ; " << std::endl ; 
    file << std::endl ; 
    file << "  out[12] = dx ;" << std::endl ; 
    file << "  out[13] = dy ;" << std::endl ; 
    file << "  out[14] = dz ;" << std::endl ; 
    file << "  out[15] = 1 ; " << std::endl ; 
    file << std::endl ; 
    file << "  return out ; " << std::endl ; 
    file << "  };" << std::endl ; 

    // rotation Y
    file << "mat44.createRotateYMatrix = function( angle_rad )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  var angle = deg_to_rad( angle_rad ) ;" << std::endl ; 
	file << "  var s = Math.sin( angle_rad ) ;" << std::endl ; 
	file << "  var c = Math.cos( angle_rad ) ; " << std::endl ; 
	file << "  var out = new Float32Array( 16 ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  out[0] = c ;" << std::endl ; 
	file << "  out[1] = 0 ;" << std::endl ; 
	file << "  out[2] = s ;" << std::endl ; 
	file << "  out[3] = 0 ;" << std::endl ; 
	file << std::endl ;
	file << "  out[4] = 0 ;" << std::endl ; 
	file << "  out[5] = 1 ;" << std::endl ; 
	file << "  out[6] = 0 ;" << std::endl ;  
	file << "  out[7] = 0 ;" << std::endl ; 
	file << std::endl ;
	file << "  out[8] = -s ;" << std::endl ;  
	file << "  out[9] = 0 ;" << std::endl ; 
	file << "  out[10] = c ;" << std::endl ; 
	file << "  out[11] = 0 ;" << std::endl ;  
	file << std::endl ;
	file << "  out[12] = 0 ;" << std::endl ; 
	file << "  out[13] = 0 ;" << std::endl ; 
	file << "  out[14] = 0 ;" << std::endl ; 
	file << "  out[15] = 1 ;" << std::endl ;  
	file << std::endl ;
	file << "  return out ;" << std::endl ;  
	file << "};" << std::endl ; 

	file.close() ; 
}

void writeMainWebGLFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "var gl ;" << std::endl;
	file << std::endl ; 
	file << "var pointShaderProgram ;" << std::endl ; 
	file << "var lineShaderProgram ;" << std::endl ; 
	file << "var surfaceShaderProgram ;" << std::endl ; 
	file << std::endl ; 
	file << "var pCloud ;" << std::endl	;
	file << std::endl ; 

	// init canvas 	
	file << "function initGL()" << std::endl ; 
	file << "{" << std::endl ;
	file << "  var canvas = document.getElementById(\"glCanvas\") ;" << std::endl ; 
	file << "  gl         = canvas.getContext(\"webgl\") || canvas.getContext(\"experimental-webgl\") ;" << std::endl ;  
	file << "  if( ! gl )" << std::endl ;
	file << "  {" << std::endl ; 
	file << "    alert(\"could not load WebGL context\");" << std::endl ; 
	file << "    return ;" << std::endl ; 
	file << "  }" << std::endl ;
	file << "  setup( ) ;" << std::endl ;  
	file << "  update( ) ;" << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ; 

	// setup rendering 
	file << "function setup( )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  gl.clearColor( 0.0 , 0.0 , 0.0 , 1.0 ) ;" << std::endl ; 
	file << "  gl.enable( gl.DEPTH_TEST ) ;" << std::endl ; 
	file << "  gl.depthFunc( gl.LEQUAL ) ;" << std::endl ; 
	file << std::endl ;
	file << "  pointShaderProgram = new Shader( gl , \"point_vShader\" , \"point_fShader\" ) ;" << std::endl ; 
	file << "  pCloud = new PointCloud( gl , pointShaderProgram , cloud_data ) ;" << std::endl ; 
	file << "  surfaceShaderProgram = new Shader( gl , \"image_vShader\" , \"image_fShader\" ) ;" << std::endl ; 
	file << "  cams   = new Array( cameras.length ) ;" << std::endl ; 
	file << "  for( var i = 0 ; i < cameras.length ; ++i )" << std::endl ; 
	file << "  {" << std::endl; 
	file << "    cams[i] = new Camera( gl , surfaceShaderProgram , pointShaderProgram , cameras[i].position , cameras[i].imagePlane , cameras[i].imageName ) ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ; 

	// update scene (moving camera lights ...)
	file << "function update()" << std::endl ; 
	file << "{" << std::endl;
	file << "  requestAnimFrame( update );" << std::endl ; 
	file << "  var tra     = mat44.createTranslationMatrix( - scene_center[0] , - scene_center[1] , - scene_center[2] ) ;" << std::endl ; 
	file << "  var inv_tra = mat44.createTranslationMatrix( scene_center[0] , scene_center[1] , scene_center[2] );" << std::endl;
	file << "  var rot     = mat44.createRotateYMatrix( 0.005 ) ;" << std::endl ; 
	file << "  var geom    = mat44.mul( tra , mat44.mul( rot , inv_tra ) );" << std::endl ; 
	file << "  cameraModelView = mat44.mul( geom , cameraModelView ) ;" << std::endl ; 
	file << "  render() ;" << std::endl; 
	file << "}" << std::endl ; 

	// render scene 
	file << "function render( )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  gl.clear( gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT ) ;" << std::endl ;
	file << "  pointShaderProgram.Enable( gl ) ;" << std::endl ; 
	file << "  var mv_loc = gl.getUniformLocation( pointShaderProgram.shad , \"uModelViewMatrix\" ) ;" << std::endl ; 
	file << "  var pj_loc = gl.getUniformLocation( pointShaderProgram.shad , \"uProjectionMatrix\" ) ;" << std::endl ; 
	file << "  gl.uniformMatrix4fv( mv_loc , false , cameraModelView ) ;" << std::endl ; 
	file << "  gl.uniformMatrix4fv( pj_loc , false , cameraProjection ) ;" << std::endl ;
	file << "  pCloud.render( gl ) ; " << std::endl ; 
	file << "  for( var i = 0 ; i < cameras.length ; ++i )" << std::endl ; 
	file << "  {" << std::endl ; 
	file << "    cams[i].renderLines( gl );" << std::endl;
	file << "  }" << std::endl; 
	file << "  surfaceShaderProgram.Enable( gl ) ; " << std::endl ;
	file << "  var mv_loc = gl.getUniformLocation( surfaceShaderProgram.shad , \"uModelViewMatrix\" ) ;" << std::endl ; 
	file << "  var pj_loc = gl.getUniformLocation( surfaceShaderProgram.shad , \"uProjectionMatrix\" ) ;" << std::endl ; 
	file << "  gl.uniformMatrix4fv( mv_loc , false , cameraModelView ) ;" << std::endl ; 
	file << "  gl.uniformMatrix4fv( pj_loc , false , cameraProjection ) ;" << std::endl ;	 
	file << "  for( var i = 0 ; i < cameras.length ; ++i )" << std::endl ; 
	file << "  {" << std::endl; 
	file << "    cams[i].render( gl ) ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << "}" << std::endl ; 



	file.close() ;
}

void WriteCameraFile( Document & doc , const std::string & file_name ) 
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	// camera class 
	file<< "function Camera( aGLContext , aShader , aShaderLine , aPosition , aPlanePos , aFileName )" << std::endl ; 
	file<< "{" << std::endl ; 
	file<< "  this.pos           = aPosition ;" << std::endl ; 
	file<< "  this.planePos      = aPlanePos ;" << std::endl ;  
	file<< "  this.shader        = aShader ;" << std::endl ;
	file<< "  this.shaderLines   = aShaderLine;" << std::endl ;  
	file<< "  this.imageFileName = aFileName ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // array = 2 triangles : 0 1 2 - 0 2 3" << std::endl ; 
	file<< "  var buff_img = new Float32Array( 30 ) ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 0 " << std::endl ; 
	file<< "  buff_img[0] = aPlanePos[0] ;" << std::endl ; 
	file<< "  buff_img[1] = aPlanePos[1] ;" << std::endl ; 
	file<< "  buff_img[2] = aPlanePos[2] ;" << std::endl ; 
	file<< "  buff_img[3] = 0.0 ;" << std::endl ; 
	file<< "  buff_img[4] = 0.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 1 " << std::endl ; 
	file<< "  buff_img[5] = aPlanePos[3] ;" << std::endl ; 
	file<< "  buff_img[6] = aPlanePos[4] ;" << std::endl ; 
	file<< "  buff_img[7] = aPlanePos[5] ;" << std::endl ; 
	file<< "  buff_img[8] = 1.0 ;" << std::endl ; 
	file<< "  buff_img[9] = 0.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 2" << std::endl ; 
	file<< "  buff_img[10] = aPlanePos[6] ;" << std::endl ; 
	file<< "  buff_img[11] = aPlanePos[7] ;" << std::endl ; 
	file<< "  buff_img[12] = aPlanePos[8] ;" << std::endl ; 
	file<< "  buff_img[13] = 1.0 ;" << std::endl ; 
	file<< "  buff_img[14] = 1.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 0 " << std::endl ; 
	file<< "  buff_img[15] = aPlanePos[0] ;" << std::endl ; 
	file<< "  buff_img[16] = aPlanePos[1] ;" << std::endl ; 
	file<< "  buff_img[17] = aPlanePos[2] ;" << std::endl ; 
	file<< "  buff_img[18] = 0.0 ;" << std::endl ; 
	file<< "  buff_img[19] = 0.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 2" << std::endl ; 
	file<< "  buff_img[20] = aPlanePos[6] ;" << std::endl ; 
	file<< "  buff_img[21] = aPlanePos[7] ;" << std::endl ; 
	file<< "  buff_img[22] = aPlanePos[8] ;" << std::endl ; 
	file<< "  buff_img[23] = 1.0 ;" << std::endl ; 
	file<< "  buff_img[24] = 1.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  // 3 " << std::endl ; 
	file<< "  buff_img[25] = aPlanePos[9] ;" << std::endl ; 
	file<< "  buff_img[26] = aPlanePos[10] ;" << std::endl ; 
	file<< "  buff_img[27] = aPlanePos[11] ;" << std::endl ; 
	file<< "  buff_img[28] = 0.0 ;" << std::endl ; 
	file<< "  buff_img[29] = 1.0 ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  var buff_lines = new Float32Array( 24 )" << std::endl ; 
	file<< "  buff_lines[0] = aPosition[0] ;" << std::endl ; 
	file<< "  buff_lines[1] = aPosition[1] ;" << std::endl ; 
	file<< "  buff_lines[2] = aPosition[2] ;" << std::endl ; 

	file<< "  buff_lines[3] = aPlanePos[0] ;" << std::endl ; 
	file<< "  buff_lines[4] = aPlanePos[1] ;" << std::endl ; 
	file<< "  buff_lines[5] = aPlanePos[2] ;" << std::endl ; 

	file<< "  buff_lines[6] = aPosition[0] ;" << std::endl ; 
	file<< "  buff_lines[7] = aPosition[1] ;" << std::endl ; 
	file<< "  buff_lines[8] = aPosition[2] ;" << std::endl ; 

	file<< "  buff_lines[9] = aPlanePos[3] ;" << std::endl ; 
	file<< "  buff_lines[10] = aPlanePos[4] ;" << std::endl ; 
	file<< "  buff_lines[11] = aPlanePos[5] ;" << std::endl ; 

	file<< "  buff_lines[12] = aPosition[0] ;" << std::endl ; 
	file<< "  buff_lines[13] = aPosition[1] ;" << std::endl ; 
	file<< "  buff_lines[14] = aPosition[2] ;" << std::endl ; 

	file<< "  buff_lines[15] = aPlanePos[6] ;" << std::endl ; 
	file<< "  buff_lines[16] = aPlanePos[7] ;" << std::endl ; 
	file<< "  buff_lines[17] = aPlanePos[8] ;" << std::endl ; 


	file<< "  buff_lines[18] = aPosition[0] ;" << std::endl ; 
	file<< "  buff_lines[19] = aPosition[1] ;" << std::endl ; 
	file<< "  buff_lines[20] = aPosition[2] ;" << std::endl ; 

	file<< "  buff_lines[21] = aPlanePos[9] ;" << std::endl ; 
	file<< "  buff_lines[22] = aPlanePos[10] ;" << std::endl ; 
	file<< "  buff_lines[23] = aPlanePos[11] ;" << std::endl ; 

	file<< "  // create shader" << std::endl ; 
	file<< "  this.vbo = aGLContext.createBuffer() ; " << std::endl ; 
	file<< "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;" << std::endl ; 
	file<< "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , buff_img , aGLContext.STATIC_DRAW ) ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  this.vboLines = aGLContext.createBuffer() ; " << std::endl ; 
	file<< "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vboLines ) ;" << std::endl;
	file<< "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , buff_lines , aGLContext.STATIC_DRAW ) ;" << std::endl ; 
	file<< "  // create texture file " << std::endl ; 
	file<< "  var tex       = gl.createTexture() ; " << std::endl ; 
	file<< "  var img = new Image() ;" << std::endl ; 
	file<< "  this.tex = tex ; " << std::endl ; 
	file<< "  this.img = img ; " << std::endl ; 
	file<< "  img.onload    = function() { handleTexture( aGLContext , img , tex ); } ;" << std::endl ;
  	file<< "  img.src       = aFileName ;" << std::endl ; 
  	file<< "  this.nbelt     = 6 ; " << std::endl ; 
  	file<< "}" << std::endl ; 
  	file<< std::endl ; 
  	file<< "function handleTexture( aGLContext , anImage , aTexture )" << std::endl ; 
  	file<< "{" << std::endl ; 
  	file<< "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, aTexture);" << std::endl ; 
	file<< "  aGLContext.texImage2D(aGLContext.TEXTURE_2D, 0, aGLContext.RGBA, aGLContext.RGBA, aGLContext.UNSIGNED_BYTE, anImage);" << std::endl ; 
	file<< "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_MAG_FILTER, aGLContext.LINEAR);" << std::endl ; 
	file<< "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_MIN_FILTER, aGLContext.LINEAR);" << std::endl ; 
	file<< "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_WRAP_S, aGLContext.CLAMP_TO_EDGE);" << std::endl ; 
	file<< "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_WRAP_T, aGLContext.CLAMP_TO_EDGE);" << std::endl ; 
//	file<< "  aGLContext.texParameteri(aGLContext.TEXTURE_2D, aGLContext.TEXTURE_MIN_FILTER, aGLContext.LINEAR_MIPMAP_NEAREST);" << std::endl ; 
//	file<< "  aGLContext.generateMipmap(aGLContext.TEXTURE_2D);" << std::endl ; 
	file<< "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, null);" << std::endl ; 
	file<< "}" << std::endl ; 
	file<< std::endl ; 
	file<< "Camera.prototype.render = function( aGLContext )" << std::endl ; 
	file<< "{" << std::endl ; 
	file<< "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  aGLContext.enableVertexAttribArray( this.shader.attribPos ) ;" << std::endl ; 
	file<< "  aGLContext.enableVertexAttribArray( this.shader.attribTex ) ;" << std::endl ; 
	file<< std::endl ; 
	file<< "  aGLContext.vertexAttribPointer( this.shader.attribPos , 3 , aGLContext.FLOAT , false , 20 , 0 ) ;" << std::endl ; 
	file<< "  aGLContext.vertexAttribPointer( this.shader.attribTex , 2 , aGLContext.FLOAT , false , 20 , 12 ) ;" << std::endl ; 
	file<< std::endl ; 
	file << "  aGLContext.activeTexture(aGLContext.TEXTURE0);" << std::endl ; 
	file << "  aGLContext.bindTexture(aGLContext.TEXTURE_2D, this.tex);" << std::endl ; 
  	file << "  aGLContext.uniform1i(aGLContext.getUniformLocation(this.shader.shad, \"uSampler\"), 0);" << std::endl ; 
  	file << std::endl ; 
	file<< "  aGLContext.drawArrays( aGLContext.TRIANGLES , 0 , this.nbelt ) ;" << std::endl ; 
	file<< "}" << std::endl ; 
	file<< std::endl ; 
	file<< "Camera.prototype.renderLines = function( aGLContext )" << std::endl ; 
	file<< "{" << std::endl ; 
	file<< "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vboLines ) ;" << std::endl ;
	file<< "  aGLContext.enableVertexAttribArray( this.shaderLines.attribPos ) ;" << std::endl ; 
	file<< "  aGLContext.vertexAttribPointer( this.shaderLines.attribPos , 3 , aGLContext.FLOAT , false , 12 , 0 ) ;" << std::endl ; 
	file<< "  aGLContext.drawArrays( aGLContext.LINES , 0 , 8 ) ;" << std::endl ; 
	file<< "}" << std::endl ;



	file.close() ; 
}

void writeShaderFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "function Shader( aGLContext , aVertexID , aFragID )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  this.shad = setupShader( aGLContext , aVertexID , aFragID ) ; " << std::endl ; 
	file << std::endl ; 
	file << "  this.attribPos = aGLContext.getAttribLocation( this.shad , \"aPosition\" ) ;" << std::endl ; 
	file << "  this.attribTex = aGLContext.getAttribLocation( this.shad , \"aTexCoord\" ) ;" << std::endl ; 
    file << "  this.attribNor = aGLContext.getAttribLocation( this.shad , \"aNormal\" ) ;" << std::endl ;  
  	file << "  this.attribCol = aGLContext.getAttribLocation( this.shad , \"aColor\" ) ;" << std::endl ; 
  	file << "}" << std::endl ; 
  	file << std::endl ; 
  	file << "Shader.prototype.Enable = function( aGLContext )" << std::endl ; 
  	file << "{" << std::endl ; 
  	file << "  aGLContext.useProgram( this.shad ) ;" << std::endl ; 
  	file << "}" << std::endl ; 
  	file << std::endl ; 
  	file << "Shader.prototype.Disable = function( aGLContext )" << std::endl ;
  	file << "{" << std::endl ; 
  	file << "  aGLContext.useProgram( 0 ) ;" << std::endl ; 
  	file << "}" << std::endl ;  
  	file << std::endl ; 
	file << "function compileShader( aGLContext , aShaderSource , aShaderType )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  var shad = aGLContext.createShader( aShaderType ) ;" << std::endl ; 
	file << "  aGLContext.shaderSource( shad , aShaderSource ) ;" << std::endl ; 
	file << "  aGLContext.compileShader( shad ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  var ok = aGLContext.getShaderParameter( shad , aGLContext.COMPILE_STATUS ) ;" << std::endl ;
	file << "  if( ! ok )" << std::endl ; 
	file << "  {" << std::endl ; 
	file << "    throw \"could not compile shader : \" + aGLContext.getShaderInfoLog( shad ) ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << std::endl ; 
	file << "  return shad ;" << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "function setupShader( aGLContext , aVertexID , aFragID )" << std::endl ; 
	file << "{" << std::endl ;
	file << "  var vElt = document.getElementById( aVertexID ) ;" << std::endl ; 
	file << "  if( ! vElt )" << std::endl ; 
	file << "  {" << std::endl ;
	file << "    throw \"could not find vertex shader\" + aVertexID ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << "  var vStr  = vElt.text ;" << std::endl ;
	file << "  var vShad = compileShader( aGLContext , vStr , aGLContext.VERTEX_SHADER ) ;" << std::endl ; 
	file << "  var fElt = document.getElementById( aFragID ) ;" << std::endl ; 
	file << "  if( ! fElt )" << std::endl ; 
	file << "  {" << std::endl ; 
	file << "    throw \"could not find fragment shader \" + aFragID ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << "  var fStr  = fElt.text ;" << std::endl ; 
	file << "  var fShad = compileShader( aGLContext , fStr , aGLContext.FRAGMENT_SHADER ) ;" << std::endl ; 
	file << "  var pgm = aGLContext.createProgram( ) ;" << std::endl ; 
	file << "  aGLContext.attachShader( pgm , vShad ) ;" << std::endl ; 
	file << "  aGLContext.attachShader( pgm , fShad ) ;" << std::endl ; 
	file << "  aGLContext.linkProgram( pgm ) ;" << std::endl ; 
	file << "  var ok = aGLContext.getProgramParameter( pgm , aGLContext.LINK_STATUS ) ;" << std::endl ; 
	file << "  if( ! ok )" << std::endl ; 
	file << "  {" << std::endl ; 
	file << "    throw \"pgm failed to link\" + aGLContext.getProgramInfoLog( pgm ) ;" << std::endl ; 
	file << "  }" << std::endl ; 
	file << "  return pgm ;" << std::endl; 
	file << "}" << std::endl ; 

	file.close() ; 
}

void writePointCloudFile( Document & doc , const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	file << "function PointCloud( aGLContext , aShader , aPointData )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  this.pointData = aPointData ;" << std::endl ; 
	file << "  this.shader    = aShader ;" << std::endl ; 
	file << "  this.nbelt     = aPointData.length / 6 ;" << std::endl ; 
	file << std::endl ; 
	file << "  this.vbo = aGLContext.createBuffer( ) ;" << std::endl ; 
	file << "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;" << std::endl ; 
	file << "  aGLContext.bufferData( aGLContext.ARRAY_BUFFER , aPointData , aGLContext.STATIC_DRAW ) ;" << std::endl ; 
	file << "}" << std::endl ; 
	file << std::endl ; 
	file << "PointCloud.prototype.render = function( aGLContext )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  aGLContext.bindBuffer( aGLContext.ARRAY_BUFFER , this.vbo ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  aGLContext.enableVertexAttribArray( this.shader.attribPos ) ;" << std::endl ; 
	file << "  aGLContext.enableVertexAttribArray( this.shader.attribCol ) ;" << std::endl ; 
	file << std::endl ;
	file << "  aGLContext.vertexAttribPointer( this.shader.attribPos , 3 , aGLContext.FLOAT , false , 24 , 0 ) ;" << std::endl ;
	file << "  aGLContext.vertexAttribPointer( this.shader.attribCol , 3 , aGLContext.FLOAT , false , 24 , 12 ) ;" << std::endl ; 
	file << std::endl ; 
	file << "  aGLContext.drawArrays( aGLContext.POINTS , 0 , this.nbelt ) ;" << std::endl ; 
	file << "}" << std::endl ; 

	file.close() ; 
}

void writeUtilFile( const std::string & file_name )
{
	std::ofstream file( file_name.c_str() ) ;
	if( ! file.good() )
	{
		std::cerr << "Error writing file : \"" << file_name << "\"" << std::endl ; 
		exit( EXIT_FAILURE ) ; 
	}

	// requestAnimFrame 
	file << "requestAnimFrame = (function()" << std::endl ; 
	file << "{" << std::endl ;
	file << "  return window.requestAnimationFrame ||" << std::endl ; 
	file << "   window.webkitRequestAnimationFrame ||" << std::endl ;
	file << "   window.mozRequestAnimationFrame ||" << std::endl ; 
	file << "   window.oRequestAnimationFrame ||" << std::endl ; 
	file << "   window.msRequestAnimationFrame ||" << std::endl ; 
	file << "   function(/* function FrameRequestCallback */ callback, /* DOMElement Element */ element) {" << std::endl ; 
	file << "   window.setTimeout(callback, 1000/60);" << std::endl ; 
	file << " };" << std::endl ; 
	file << "})();" << std::endl ; 

	// deg_to_rad
	file << "function deg_to_rad( arad )" << std::endl ; 
	file << "{" << std::endl ; 
	file << "  return ((arad) * 0.017453292519943295769236907684886127134428718885417254560971914401710091146034494436822415696345094823 ) ;" << std::endl ; 
	file << "}" << std::endl ; 
	
	file.close() ; 
}

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


void exportProjectToWebGL( Document & doc , const std::string & SfM_outputDir , const std::string & outFolder )
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

	// write shader js file
	const std::string shader_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "shader.js" ; 
	writeShaderFile( shader_js_file_name ) ; 

	// write camera webgl file
	const std::string camera_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "cameras.js" ; 
	WriteCameraFile( doc , camera_js_file_name ) ; 

	// write pointcloud webgl file 
	const std::string pcloud_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "cloud.js" ; 
	writePointCloudFile( doc , pcloud_js_file_name ) ; 

	const std::string scene_data_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "sceneData.js" ; 
	writeSceneDataFile( doc , scene_data_file_name ) ; 

	// write main webgl file 
	const std::string main_js_file_name = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "scripts" ) + "main.js" ; 
	writeMainWebGLFile( main_js_file_name ) ; 

	// write images files
	const std::string image_folder = stlplus::folder_append_separator(stlplus::folder_append_separator( outFolder ) + "images" );
	writeImagesFiles( doc , SfM_outputDir , image_folder ) ;

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


	
	return EXIT_SUCCESS ; 
}