// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Main webgl object : contains the gl context, the shaders and the current view and projection matrices */
RenderContext = function( canvas )
{
  this.m_canvas = canvas ; 

  this.initGLContext() ;
  this.initShaders() ; 
}

/* Initialize the webgl context */
RenderContext.prototype.initGLContext = function( )
{
  try
  {
    this.m_gl = this.m_canvas.getContext( "webgl") || this.m_canvas.getContext( "experimental-webgl" ) ;
  }
  catch(e)
  {

  }

  if( ! this.m_gl )
  {
    alert( "Error Initialization of the gl context" ) ; 
  }
  this.setClearColor( 0.0 , 0.0 , 0.0 , 1.0 ) ; 
}

/* Initialize the shaders */
RenderContext.prototype.initShaders = function()
{
  // Load shaders 
  var trackballVertexShaderElemID    = document.getElementById( "line_v_shad" ) ;
  var trackballVertexShaderContent   = trackballVertexShaderElemID.text ; 

  var trackballFragmentShaderElemID  = document.getElementById( "line_f_shad" ) ;
  var trackballFragmentShaderContent = trackballFragmentShaderElemID.text ;

  this.m_trackball_shad = new Shader( trackballVertexShaderContent , trackballFragmentShaderContent , this ) ; 


  var pointVertexShaderElemID = document.getElementById( "point_v_shad" ) ;
  var pointVertexShaderContent = pointVertexShaderElemID.text ;

  var pointFragmentShaderElemID = document.getElementById( "point_f_shad" ) ;
  var pointFragmentShaderContent = pointFragmentShaderElemID.text ; 

  this.m_point_shad = new Shader( pointVertexShaderContent , pointFragmentShaderContent , this ) ; 
}

/* Set clear color */
RenderContext.prototype.setClearColor = function( aR , aG , aB , aA )
{
  this.m_gl.clearColor( aR , aG , aB , aA ) ; 
}

/* Get the webgl context */
RenderContext.prototype.getGLContext = function()
{
  return this.m_gl ;
}

/* Get the trackball shader */
RenderContext.prototype.getTrackballShader = function()
{
  return this.m_trackball_shad ;
}

RenderContext.prototype.getPointcloudShader = function()
{
  return this.m_point_shad ; 
}

/* Set the current view matrix */
RenderContext.prototype.setCurrentViewMatrix = function( aViewMat )
{
  this.m_view_mat = aViewMat ;
}

/* Set the current projection matrix */
RenderContext.prototype.setCurrentProjectionMatrix = function( aProjMat )
{
  this.m_proj_mat = aProjMat ; 
}

/* Get the current view matrix */
RenderContext.prototype.getCurrentViewMatrix = function( )
{
  return this.m_view_mat ;
}

/* Get the current projection matrix */
RenderContext.prototype.getCurrentProjectionMatrix = function()
{
  return this.m_proj_mat ; 
}
