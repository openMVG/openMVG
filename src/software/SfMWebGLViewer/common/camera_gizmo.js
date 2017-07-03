/* Create a gizmo */
CameraGizmo = function( aPosition , anImagePlane , anImage , aRenderContext )
{
  this.m_position    = aPosition;
  this.m_show        = true;
  this.m_image_plane = anImagePlane;
  this.m_image       = anImage;
  this.m_nb_vert     = 16;
  this.m_orient = new DualQuaternion();
  this.m_scale = 1.0;

  this.initGLData( aRenderContext );
}

/* Set visibility of the trackball */
CameraGizmo.prototype.setVisible = function( aVal )
{
  this.m_show = aVal;
}

/* Indicate if the trackball is visible */
CameraGizmo.prototype.isVisible = function( )
{
  return this.m_show;
}

/* Rotate a the gizmo */
CameraGizmo.prototype.rotate = function( aQuat , aCenter )
{
  var inv = Vector.negate( aCenter );

  var dqv = new DualQuaternion();
  dqv.setFromTranslationVector( aCenter );
  var dqinv = new DualQuaternion();
  dqinv.setFromTranslationVector( inv );
  var dqq = new DualQuaternion();
  dqq.setFromRotationQuaternion( aQuat );

  this.m_orient = dqv.mul( dqq.mul( dqinv.mul( this.m_orient ) ) );
}

/* Translate the gizmo */
CameraGizmo.prototype.translate = function( aVector )
{
  var dqv = new DualQuaternion();
  dqv.setFromTranslationVector( aVector );

  this.m_orient = dqv.mul( this.m_orient );
}

/* Reset to it's initial position */
CameraGizmo.prototype.reset = function()
{
  this.m_orient = new DualQuaternion();
  this.m_scale = 1.0;
}

/* Set scaling factor */
CameraGizmo.prototype.setScale = function( aScaleValue )
{
  this.m_scale = aScaleValue;
}

/* Draw the gizmo */
CameraGizmo.prototype.initGLData = function( aRenderContext )
{
  var gl = aRenderContext.getGLContext();

  var lineData = new Common.ARRAY_TYPE( 2 * 6 * 8 ); // 2 points per line; 6 coord per point; 8 lines

  // From points to imagePlane
  for (var lineId = 0; lineId < 4; ++lineId )
  {
    // Position of first point
    lineData[ 12 * lineId ]     = this.m_position[ 0 ];
    lineData[ 12 * lineId + 1 ] = this.m_position[ 1 ];
    lineData[ 12 * lineId + 2 ] = this.m_position[ 2 ];
    // Color of first point (white)
    lineData[ 12 * lineId + 3 ] = 1.0;
    lineData[ 12 * lineId + 4 ] = 1.0;
    lineData[ 12 * lineId + 5 ] = 1.0;

    // Position of second point
    lineData[ 12 * lineId + 6 ] = this.m_image_plane[ 3 * lineId ];
    lineData[ 12 * lineId + 7 ] = this.m_image_plane[ 3 * lineId + 1 ];
    lineData[ 12 * lineId + 8 ] = this.m_image_plane[ 3 * lineId + 2 ];
    // Color of second point
    lineData[ 12 * lineId + 9 ]  = 1.0;
    lineData[ 12 * lineId + 10 ] = 1.0;
    lineData[ 12 * lineId + 11 ] = 1.0;
  }

  var pad = 12 * 4;

  // Between images planes
  for (lineId = 0; lineId < 4; ++lineId )
  {
    var pt_id      = lineId;
    var pt_id_next = (lineId + 1) % 4;

    // Position of first point
    lineData[ pad + lineId * 12 ]     = this.m_image_plane[ 3 * pt_id ];
    lineData[ pad + lineId * 12 + 1 ] = this.m_image_plane[ 3 * pt_id + 1 ];
    lineData[ pad + lineId * 12 + 2 ] = this.m_image_plane[ 3 * pt_id + 2 ];
    // Color of first point
    lineData[ pad + lineId * 12 + 3 ] = 1.0;
    lineData[ pad + lineId * 12 + 4 ] = 1.0;
    lineData[ pad + lineId * 12 + 5 ] = 1.0;

    // Position of second point
    lineData[ pad + lineId * 12 + 6 ] = this.m_image_plane[ 3 * pt_id_next ];
    lineData[ pad + lineId * 12 + 7 ] = this.m_image_plane[ 3 * pt_id_next + 1 ];
    lineData[ pad + lineId * 12 + 8 ] = this.m_image_plane[ 3 * pt_id_next + 2 ];
    // Color of second point
    lineData[ pad + lineId * 12 + 9 ]  = 1.0;
    lineData[ pad + lineId * 12 + 10 ] = 1.0;
    lineData[ pad + lineId * 12 + 11 ] = 1.0;
  }

  this.m_vbo = gl.createBuffer();
  gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo );
  gl.bufferData( gl.ARRAY_BUFFER , lineData , gl.STATIC_DRAW );
}

/* Draw the camera gizmo */
CameraGizmo.prototype.draw = function( aRenderContext )
{
  if (this.m_show )
  {
    var shad = aRenderContext.getTrackballShader();
    var gl = aRenderContext.getGLContext();

    shad.enable();
    // Set common matrices
    var dqv = new DualQuaternion();
    dqv.setFromTranslationVector( Vector.negate( this.m_position ) );
    var scale = Matrix4.createUniformScale( this.m_scale );
    var dqvinv = new DualQuaternion();
    dqvinv.setFromTranslationVector( this.m_position );


    var model_mat = Matrix4.transpose( Matrix4.mul( this.m_orient.toMatrix() , Matrix4.mul( dqvinv.toMatrix() , Matrix4.mul( scale , dqv.toMatrix() ) ) ) );
    var mvp = Matrix4.mul( model_mat , Matrix4.mul( aRenderContext.getCurrentViewMatrix() , aRenderContext.getCurrentProjectionMatrix() ) );
  shad.setModelViewProjectionMatrix( mvp );
;
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
