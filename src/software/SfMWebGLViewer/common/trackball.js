// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Class for manipulation and drawing of a trackball */
Trackball = function( trackballResolution , aRenderContext ) 
{
  this.m_resolution = trackballResolution;
  this.m_show       = true; 
  this.m_nb_vert    = 2 * 3 * this.m_resolution; 

  this.m_orient     = new DualQuaternion(); 
  this.m_radius     = 1.0; 

  this.initGLData( aRenderContext );
}

/* Create vertex buffer data */
Trackball.prototype.initGLData = function( aRenderContext )
{
  var trackBallNbComponentPerVert = 6; // 3 for color + 3 for position 
  var trackBallData = new Common.ARRAY_TYPE( this.m_nb_vert * trackBallNbComponentPerVert );

  // X - data 
  var pad = 0; 
  for (var i = 0; i < this.m_resolution; ++i )
  {
    var theta = i * 2.0 * Math.PI / this.m_resolution;
    var theta_ip = ((i+1) % this.m_resolution ) * 2.0 * Math.PI / this.m_resolution;

    var cs = Math.cos( theta );
    var ss = Math.sin( theta );

    var cs_ip = Math.cos( theta_ip );
    var ss_ip = Math.sin( theta_ip );

    // first vertex 
    trackBallData[ 12 * i ] = 0.0;
    trackBallData[ 12 * i + 1 ] = cs;
    trackBallData[ 12 * i + 2 ] = ss;
    // - RGB 
    trackBallData[ 12 * i + 3 ] = 1.0;
    trackBallData[ 12 * i + 4 ] = 0.0;
    trackBallData[ 12 * i + 5 ] = 0.0; 

    // second vertex 
    trackBallData[ 12 * i + 6 ] = 0.0;
    trackBallData[ 12 * i + 7 ] = cs_ip;
    trackBallData[ 12 * i + 8 ] = ss_ip;
    // RGB 
    trackBallData[ 12 * i + 9 ]  = 1.0;
    trackBallData[ 12 * i + 10 ] = 0.0;
    trackBallData[ 12 * i + 11 ] = 0.0; 
  }
  pad += 12 * this.m_resolution; 

  // Y - data 
  for (var i = 0; i < this.m_resolution; ++i )
  {
    var theta = i * 2.0 * Math.PI / this.m_resolution;
    var theta_ip = ((i+1) % this.m_resolution ) * 2.0 * Math.PI / this.m_resolution;

    var cs = Math.cos( theta );
    var ss = Math.sin( theta );

    var cs_ip = Math.cos( theta_ip );
    var ss_ip = Math.sin( theta_ip );

    // first vertex 
    trackBallData[ pad + 12 * i ] = cs;
    trackBallData[ pad + 12 * i + 1 ] = 0.0;
    trackBallData[ pad + 12 * i + 2 ] = ss;
    // - RGB 
    trackBallData[ pad + 12 * i + 3 ] = 0.0;
    trackBallData[ pad + 12 * i + 4 ] = 1.0;
    trackBallData[ pad + 12 * i + 5 ] = 0.0; 

    // second vertex 
    trackBallData[ pad + 12 * i + 6 ] = cs_ip;
    trackBallData[ pad + 12 * i + 7 ] = 0.0;
    trackBallData[ pad + 12 * i + 8 ] = ss_ip;
    // RGB 
    trackBallData[ pad + 12 * i + 9 ]  = 0.0;
    trackBallData[ pad + 12 * i + 10 ] = 1.0;
    trackBallData[ pad + 12 * i + 11 ] = 0.0; 
  }
  pad += 12 * this.m_resolution;
  // Z - data 
  for (var i = 0; i < this.m_resolution; ++i )
  {
    var theta = i * 2.0 * Math.PI / this.m_resolution;
    var theta_ip = ((i+1) % this.m_resolution ) * 2.0 * Math.PI / this.m_resolution;

    var cs = Math.cos( theta );
    var ss = Math.sin( theta );

    var cs_ip = Math.cos( theta_ip );
    var ss_ip = Math.sin( theta_ip );

    // first vertex 
    trackBallData[ pad + 12 * i ] = cs;
    trackBallData[ pad + 12 * i + 1 ] = ss;
    trackBallData[ pad + 12 * i + 2 ] = 0.0;
    // - RGB 
    trackBallData[ pad + 12 * i + 3 ] = 0.0;
    trackBallData[ pad + 12 * i + 4 ] = 0.0;
    trackBallData[ pad + 12 * i + 5 ] = 1.0; 

    // second vertex 
    trackBallData[ pad + 12 * i + 6 ] = cs_ip;
    trackBallData[ pad + 12 * i + 7 ] = ss_ip; 
    trackBallData[ pad + 12 * i + 8 ] = 0.0;
    // RGB 
    trackBallData[ pad + 12 * i + 9 ]  = 0.0;
    trackBallData[ pad + 12 * i + 10 ] = 0.0;
    trackBallData[ pad + 12 * i + 11 ] = 1.0; 
  }

  // Build VBO 
  var gl = aRenderContext.getGLContext(); 
  this.m_vbo =  gl.createBuffer( );
  gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo ); 
  gl.bufferData( gl.ARRAY_BUFFER , trackBallData , gl.STATIC_DRAW );
}

/* Set visibility of the trackball */
Trackball.prototype.setVisible = function( aVal )
{
  this.m_show = aVal;
}

/* Indicate if the trackball is visible */
Trackball.prototype.isVisible = function( )
{
  return this.m_show; 
}

/* Get orientation (dual) quaternion */
Trackball.prototype.getOrientation = function( )
{
  return this.m_orient;
}

/* Get set orientation (dual quaternion) */
Trackball.prototype.setOrientation = function( aNewOrientation )
{
  this.m_orient.copy( aNewOrientation );
}

/* Draw on a webgl context */
Trackball.prototype.draw = function( aRenderContext )
{ 
  if (this.m_show )
  {
    var shad = aRenderContext.getTrackballShader();
    var gl = aRenderContext.getGLContext();

    shad.enable(); 
    // Set common matrices 
    var model_mat = Matrix4.transpose( Matrix4.mul( this.m_orient.toMatrix() , Matrix4.createUniformScale( this.m_radius ) ) );
    var mvp = Matrix4.mul( Matrix4.mul( model_mat , aRenderContext.getCurrentViewMatrix() ) , aRenderContext.getCurrentProjectionMatrix() );
    shad.setModelViewProjectionMatrix( mvp ); 

    // Now render 
    gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo );

    var posLoc = shad.getPositionAttributeLocation();
    var colLoc = shad.getColorAttributeLocation();
  
    gl.enableVertexAttribArray( posLoc );
    gl.enableVertexAttribArray( colLoc ); 

    gl.vertexAttribPointer( posLoc , 3 , gl.FLOAT , false , ( ( 3 + 3 ) * 4 ) , 0 );
    gl.vertexAttribPointer( colLoc , 3 , gl.FLOAT , false , ( ( 3 + 3 ) * 4 ) , 3 * 4 ); 

    gl.drawArrays( gl.LINES , 0 , this.m_nb_vert ); 
  }
}

/* Set position of the trackball */
Trackball.prototype.setPosition = function( aPos )
{
  this.m_orient.setFromTranslationVector( aPos ); 
}

/* Set radius of the trackball */
Trackball.prototype.setRadius = function( aRad )
{
  this.m_radius = aRad; 
}

/* Get current radius */
Trackball.prototype.getRadius = function( )
{
  return this.m_radius; 
}

/* Reset orientation */
Trackball.prototype.reset = function()
{
  this.m_orient = new DualQuaternion(); 
}