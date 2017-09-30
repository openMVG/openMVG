// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Main webgl object : contains the gl context, the shaders and the current view and projection matrices */
RenderContext = function( canvas )
{
  this.m_canvas = canvas;
  this.m_current_cam = undefined;

  this.initGLContext();
  this.initShaders();

  this.m_backface_culling_active = true;
  this.m_view_dependent_point_size_active = true;
}

// Add your prefix here.
var browserPrefixes = [
  "",
  "MOZ_",
  "OP_",
  "WEBKIT_"
];

var getExtensionWithKnownPrefixes = function(gl, name) {
  for (var ii = 0; ii < browserPrefixes.length; ++ii) {
    var prefixedName = browserPrefixes[ii] + name;
    var ext = gl.getExtension(prefixedName);
    if (ext) {
      return ext;
    }
  }
};

/* Initialize the webgl context */
RenderContext.prototype.initGLContext = function( )
{
  console.log("coucou");
  try
  {
    this.m_gl = this.m_canvas.getContext( "webgl") || this.m_canvas.getContext( "experimental-webgl" );
  }
  catch (e)
  {

  }

  if( ! this.m_gl )
  {
    alert( "Error Initialization of the gl context" );
  }
  var ext = getExtensionWithKnownPrefixes( this.m_gl , "EXT_frag_depth" );

  this.setClearColor( 0.0 , 0.0 , 0.0 , 1.0 );
  this.m_gl.enable(this.m_gl.DEPTH_TEST);
}

/* Initialize the shaders */
RenderContext.prototype.initShaders = function()
{
  // Load shaders
  var trackballVertexShaderElemID    = document.getElementById( "line_v_shad" );
  var trackballVertexShaderContent   = trackballVertexShaderElemID.text;

  var trackballFragmentShaderElemID  = document.getElementById( "line_f_shad" );
  var trackballFragmentShaderContent = trackballFragmentShaderElemID.text;

  this.m_trackball_shad = new Shader( trackballVertexShaderContent , trackballFragmentShaderContent , this );


  var pointVertexShaderElemID = document.getElementById( "point_v_shad" );
  var pointVertexShaderContent = pointVertexShaderElemID.text;

  var pointFragmentShaderElemID = document.getElementById( "point_f_shad" );
  var pointFragmentShaderContent = pointFragmentShaderElemID.text;

  this.m_point_shad = new Shader( pointVertexShaderContent , pointFragmentShaderContent , this );
}

/* Set clear color */
RenderContext.prototype.setClearColor = function( aR , aG , aB , aA )
{
  this.m_gl.clearColor( aR , aG , aB , aA );
}

/* Get the webgl context */
RenderContext.prototype.getGLContext = function()
{
  return this.m_gl;
}

/* Get the trackball shader */
RenderContext.prototype.getTrackballShader = function()
{
  return this.m_trackball_shad;
}

RenderContext.prototype.getPointcloudShader = function()
{
  return this.m_point_shad;
}

/* Set the current view matrix */
RenderContext.prototype.setCurrentViewMatrix = function( aViewMat )
{
  this.m_view_mat = aViewMat;
}

/* Set the current projection matrix */
RenderContext.prototype.setCurrentProjectionMatrix = function( aProjMat )
{
  this.m_proj_mat = aProjMat;
}

/* Get the current view matrix */
RenderContext.prototype.getCurrentViewMatrix = function( )
{
  return this.m_view_mat;
}

/* Get the current projection matrix */
RenderContext.prototype.getCurrentProjectionMatrix = function()
{
  return this.m_proj_mat;
}

/**
 * @brief Set current camera
 * @param aCam a new camera
 */
RenderContext.prototype.setCurrentCamera = function( aCam )
{
  this.m_current_cam = aCam;
}

/**
 * @brief Get current camera
 * @return Current camera
 */
RenderContext.prototype.getCurrentCamera = function()
{
  return this.m_current_cam;
}


/**
 * Indicate if backface culling is active
 * @retval true if backface culling is active
 * @retval false if not active
 */
RenderContext.prototype.isBackfaceCullingActive = function()
{
  return this.m_backface_culling_active;
}

/**
 * @brief Get current state of the backface culling algorithm
 * @retval true If BC is active
 * @retval false If BC is inactive
 */
RenderContext.prototype.backfaceCullingState = function()
{
  return this.m_backface_culling_active;
}

/**
 * @brief Set backface culling state
 * @param aValue new state of the backface culling algorithm
 */
RenderContext.prototype.setBackFaceCulling = function( aValue )
{
  this.m_backface_culling_active = aValue;
}

/**
 * @brief get an array of two elements
 */
RenderContext.prototype.getScreenSize = function( aValue )
{
  var res = new Common.INT_ARRAY_TYPE(2);
  res[0] = this.m_current_cam.m_width;
  res[1] = this.m_current_cam.m_height;
  return res;
}

RenderContext.prototype.isViewDependentPointSizeActive = function( aValue )
{
  return this.m_view_dependent_point_size_active;
}

RenderContext.prototype.setViewDependentPointSize = function( aValue )
{
  this.m_view_dependent_point_size_active = aValue;
}
