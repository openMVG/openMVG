PointCloud = function( aPointPosition , aPointNormal , aPointColor , aRenderContext )
{
  this.initGLData( aPointPosition , aPointNormal , aPointColor , aRenderContext ) ;
}

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

PointCloud.prototype.draw = function( aRenderContext )
{
  var gl = aRenderContext.getGLContext() ;
  var shad = aRenderContext.getPointcloudShader() ;
  shad.enable() ; 

  gl.bindBuffer( gl.ARRAY_BUFFER , this.m_vbo ) ; 

  var mvp = Matrix4.mul( aRenderContext.getCurrentViewMatrix() , aRenderContext.getCurrentProjectionMatrix() ) ;
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
