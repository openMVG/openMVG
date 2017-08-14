// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Helper class for manipulation of a shader */
Shader = function( aVertexShaderContent , aFragmentShaderContent , aRenderContext )
{
  this.m_vert_content = aVertexShaderContent;
  this.m_frag_content = aFragmentShaderContent;
  this.m_gl = aRenderContext.getGLContext(); 

  this.initShader( this.m_vert_content , this.m_frag_content );
  this.getStandardLocations( ); 
}

/* shader compilation (vertex of fragment)*/
function compileShader( aGLContext , aShaderSource , aShaderType )
{
  var shad = aGLContext.createShader( aShaderType );
  aGLContext.shaderSource( shad , aShaderSource );
  aGLContext.compileShader( shad );

  var ok = aGLContext.getShaderParameter( shad , aGLContext.COMPILE_STATUS );
  if( ! ok )
  {
    throw "could not compile shader : " + aGLContext.getShaderInfoLog( shad );
  }

  return shad;
}

/* Compile vertex and fragment shaders then assemble in a program */
Shader.prototype.initShader = function( aVertexShaderContent , aFragmentShaderContent )
{
  this.m_vert_shad = compileShader( this.m_gl , aVertexShaderContent , this.m_gl.VERTEX_SHADER );
  this.m_frag_shad = compileShader( this.m_gl , aFragmentShaderContent , this.m_gl.FRAGMENT_SHADER ); 
  this.m_pgm = this.m_gl.createProgram( );
  this.m_gl.attachShader( this.m_pgm , this.m_vert_shad );
  this.m_gl.attachShader( this.m_pgm , this.m_frag_shad );
  this.m_gl.linkProgram( this.m_pgm ); 

  var ok = this.m_gl.getProgramParameter( this.m_pgm , this.m_gl.LINK_STATUS );
  if( ! ok )
  {
    console.log( "pgm failed to link" + this.m_gl.getProgramInfoLog( this.m_pgm ) );
  }
}

/* Activate this shader */
Shader.prototype.enable = function( )
{
  this.m_gl.useProgram( this.m_pgm );
}

/* deactivate this shader */
Shader.prototype.disable = function( )
{
  this.m_gl.useProgram( 0 ); 
}

/* Retreive standard locations names (inPos, inNor, inCol, uModelViewProjMat,uNormalMat) */
Shader.prototype.getStandardLocations = function(  )
{
  // Standard attributes 
  this.m_pos_loc = this.m_gl.getAttribLocation( this.m_pgm , "inPos" );
  this.m_nor_loc = this.m_gl.getAttribLocation( this.m_pgm , "inNor" );
  this.m_col_loc = this.m_gl.getAttribLocation( this.m_pgm , "inCol" ); 

  // Standard uniforms values 
  this.m_unif_mvp_mat    = this.m_gl.getUniformLocation( this.m_pgm , "uModelViewProjMat" );
  this.m_unif_mv_mat     = this.m_gl.getUniformLocation( this.m_pgm , "uModelViewMat" ); 
  this.m_unif_proj_mat   = this.m_gl.getUniformLocation( this.m_pgm , "uProjMat" ); 
  this.m_unif_normal_mat = this.m_gl.getUniformLocation( this.m_pgm , "uNormalMat" );
  
  this.m_unif_camera_pos = this.m_gl.getUniformLocation( this.m_pgm , "uCamPos" ); 
  this.m_unif_light_pos  = this.m_gl.getUniformLocation( this.m_pgm , "uLightPos" ); 
}

/* Indicate if the shader contains a position attribute */
Shader.prototype.hasPositionAttribute = function( )
{
  return this.m_pos_loc != -1;
}

/* Get position attribute location */
Shader.prototype.getPositionAttributeLocation = function( )
{
  return this.m_pos_loc; 
}

/* Indicate if the shader contains a normal attribute */
Shader.prototype.hasNormalAttribute = function( )
{
  return this.m_nor_loc != -1;
}

/* Get normal attribute location */
Shader.prototype.getNormalAttributeLocation = function()
{
  return this.m_nor_loc; 
}

/* Indicate if the shader contains a color attribute */
Shader.prototype.hasColorAttribute = function()
{
  return this.m_col_loc != -1; 
}

/* Get color attribute location */
Shader.prototype.getColorAttributeLocation = function()
{
  return this.m_col_loc; 
}

/* Indicate if the shader contains light position uniform */
Shader.prototype.hasLightPositionUniform = function()
{
  return this.m_unif_light_pos != -1; 
}

/* Get location of light position uniform */
Shader.prototype.getLightPositionUniform = function()
{
  return this.m_unif_light_pos; 
}

/* Indicate if the shader contains camera position uniform */
Shader.prototype.hasCameraPositionUniform = function()
{
  return this.m_unif_camera_pos != undefined; 
}

/* Get location of camera position uniform */
Shader.prototype.getCameraPositionUniform = function()
{
  return this.m_unif_camera_pos;
}

/**
 * @brief Set Camera position uniform 
 * @param aPosition The position to pass to the shader (array of 3 floats)
 */
Shader.prototype.setCameraPosition = function( aPosition )
{
  if (this.m_unif_camera_pos != undefined )
  {
    this.m_gl.uniform3f( this.m_unif_camera_pos , aPosition[0] , aPosition[1] , aPosition[2] );
  }
}

/* Get a location by it's name */
Shader.prototype.getAttribLocation = function( aName ) 
{
  return this.m_gl.getAttribLocation( this.m_pgm , aName );
}

/* Set the mvp projection matrix */
Shader.prototype.setModelViewProjectionMatrix = function( aMatrix )
{
  this.m_gl.uniformMatrix4fv( this.m_unif_mvp_mat , false , aMatrix ); 
}

/**
 * @brief Set the model view matrix 
 * @param aMatrix The new model view matrix 
 * @note this function assume there's a model view matrix uniform and no validity test is made
 */
Shader.prototype.setModelViewMatrix = function( aMatrix )
{
  this.m_gl.uniformMatrix4fv( this.m_unif_mv_mat , false , aMatrix ); 
}

/**
 * @brief Set the normal matrix to the shader 
 * @param aMatrix the new normal matrix 
 * @note this function assume there's a normal matrix uniform and no validity test is made
 */
Shader.prototype.setNormalMatrix = function( aMatrix )
{
  this.m_gl.uniformMatrix4fv( this.m_unif_normal_mat , false , aMatrix );
}

/**
 * @brief Set the projection matrix to the shader  
 * @param aMatrix the new projection matrix 
 */
Shader.prototype.setProjectionMatrix = function( aMatrix )
{
  this.m_gl.uniformMatrix4fv( this.m_unif_proj_mat , false , aMatrix ); 
}

/**
 * @brief Pass a float value to the shader 
 * @param aName Name of the uniform 
 * @param aValue Value to pass 
 * @note a test is made to validate if the uniform is present in the shader
 * @note If the uniform is not present, the function do nothing. 
 */
Shader.prototype.setFloatUniform = function( aName , aValue )
{
  var loc = this.m_gl.getUniformLocation( this.m_pgm , aName );
  if (loc != null )
  {
    this.m_gl.uniform1f( loc , aValue ); 
  }
}

/**
 * @brief Pass a bool value to the shader 
 * @param aName Name of the uniform 
 * @param aValue Value to pass 
 * @note a test is made to validate if the uniform is present in the shader
 * @note If the uniform is not present, the function do nothing. 
 */
Shader.prototype.setBoolUniform = function( aName , aValue )
{
  var loc = this.m_gl.getUniformLocation( this.m_pgm , aName );
  if (loc != null )
  {
    this.m_gl.uniform1i( loc , aValue ? 1 : 0 ); 
  }
}

/**
 * @brief Pass an integer array to a uniform 
 * @param aName Name of the uniform 
 * @param aValue Value to pass 
 * @note a test is made to validate if the uniform is present in the shader
 * @note If the uniform is not present, the function do nothing. 
 */
Shader.prototype.setIntegerArrayUniform = function( aName , aValue )
{
  var loc = this.m_gl.getUniformLocation( this.m_pgm , aName );
  if (loc != null )
  {
    this.m_gl.uniform1iv( loc , aValue );
  }
}