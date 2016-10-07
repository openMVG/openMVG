/* Class managing a point cloud */
PointCloud = function( aPointPosition , aPointNormal , aPointColor , aRenderContext )
{
  this.initGLData( aPointPosition , aPointNormal , aPointColor , aRenderContext ) ;
  this.m_orient = new DualQuaternion() ; 
  this.m_show   = true ; 
}

/* Create vbo */
PointCloud.prototype.initGLData = function( aPointPosition , aPointNormal , aPointColor , aRenderContext )
{
  this.m_nb_point = aPointPosition.length / 3 ;

  var gl = aRenderContext.getGLContext() ; 

  var tmpData = new Common.ARRAY_TYPE( this.m_nb_point * 3 * 3 ) ; // x3 -> three arrays ; x3 -> three components per array 

  for( var i = 0  ; i < this.m_nb_point  ; ++i )
  {
    tmpData[ 9 * i ]     = aPointPosition[ 3 * i ] ;
    tmpData[ 9 * i + 1 ] = aPointPosition[ 3 * i + 1 ] ;  
    tmpData[ 9 * i + 2 ] = aPointPosition[ 3 * i + 2 ] ;
    
    if( aPointNormal != undefined )
    {
      tmpData[ 9 * i + 3 ] = aPointNormal[ 3 * i ] ;
      tmpData[ 9 * i + 4 ] = aPointNormal[ 3 * i + 1 ] ;
      tmpData[ 9 * i + 5 ] = aPointNormal[ 3 * i + 2 ] ;
    }
    else 
    {
      tmpData[ 9 * i + 3 ] = 0.0 ; 
      tmpData[ 9 * i + 4 ] = 0.0 ; 
      tmpData[ 9 * i + 5 ] = 0.0 ;
    }

    if( aPointColor != undefined )
    {
      tmpData[ 9 * i + 6 ] = aPointColor[ 3 * i ] / 255.0 ;
      tmpData[ 9 * i + 7 ] = aPointColor[ 3 * i + 1 ] / 255.0 ; 
      tmpData[ 9 * i + 8 ] = aPointColor[ 3 * i + 2 ] / 255.0 ; 
    }
    else 
    {
      tmpData[ 9 * i + 6 ] = 1.0 ; 
      tmpData[ 9 * i + 7 ] = 1.0 ; 
      tmpData[ 9 * i + 8 ] = 1.0 ; 
    }
  }

  this.m_vbo = gl.createBuffer() ; 
  gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo ) ; 
  gl.bufferData( gl.ARRAY_BUFFER , tmpData , gl.STATIC_DRAW ) ;
}

/* Draw a point cloud */
PointCloud.prototype.draw = function( aRenderContext )
{
  if( this.m_show )
  {
    var gl = aRenderContext.getGLContext() ;
    var shad = aRenderContext.getPointcloudShader() ;
    shad.enable() ; 

    gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo ) ; 

    var mvp = Matrix4.mul( Matrix4.transpose(this.m_orient.toMatrix()) , Matrix4.mul( aRenderContext.getCurrentViewMatrix() , aRenderContext.getCurrentProjectionMatrix() ) ) ;
  shad.setModelViewProjectionMatrix( mvp ) ;

    var posLoc = shad.getPositionAttributeLocation() ;
    if( shad.hasNormalAttribute() )
    {
      var norLoc = shad.getNormalAttributeLocation() ;
      gl.enableVertexAttribArray( norLoc ) ;
      gl.vertexAttribPointer( norLoc , 3 , gl.FLOAT , false , ( ( 3 + 3 + 3 ) * 4 ) , ( 3 * 4 ) ) ;
    }
    var colLoc = shad.getColorAttributeLocation() ;

    gl.enableVertexAttribArray( posLoc ) ;
    gl.enableVertexAttribArray( colLoc ) ;

    gl.vertexAttribPointer( posLoc , 3 , gl.FLOAT , false , ( ( 3 + 3 + 3 ) * 4 ) , 0 ) ;
    gl.vertexAttribPointer( colLoc , 3 , gl.FLOAT , false , ( ( 3 + 3 + 3 ) * 4 ) , ( 3 + 3 ) * 4 ) ; 

    gl.drawArrays( gl.POINTS , 0 , this.m_nb_point ) ; 
  }
}

/* Rotate a point cloud */
PointCloud.prototype.rotate = function( aQuat , aCenter )
{
  var inv = Vector.negate( aCenter ) ;

  var dqv = new DualQuaternion() ;
  dqv.setFromTranslationVector( aCenter ) ;
  var dqinv = new DualQuaternion() ;
  dqinv.setFromTranslationVector( inv ) ; 
  var dqq = new DualQuaternion() ;
  dqq.setFromRotationQuaternion( aQuat ) ;

  this.m_orient = dqv.mul( dqq.mul( dqinv.mul( this.m_orient ) ) )  ;
}

/* Translate a point cloud */
PointCloud.prototype.translate = function( aVector )
{
  var dqv = new DualQuaternion() ;
  dqv.setFromTranslationVector( aVector ) ; 

  this.m_orient = dqv.mul( this.m_orient ) ; 
}

/* Set visibility of the point cloud */
PointCloud.prototype.setVisible = function( aVal )
{
  this.m_show = aVal ;
}

/* Indicate if the point cloud is visible */
PointCloud.prototype.isVisible = function( )
{
  return this.m_show ; 
}

/* Reset orientation */
PointCloud.prototype.reset = function()
{
  this.m_orient = new DualQuaternion() ; 
}