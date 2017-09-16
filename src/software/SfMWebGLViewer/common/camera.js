// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


/*
* Create a Perspective matrix 
*/
PerspectiveCamera = function( aPosition , aDestination , aUp , aFOV , aNearClip , aFarClip , aSensorWidth , aSensorHeight )
{
  this.m_pos    = aPosition;
  this.m_dir    = aDestination;
  this.m_up     = aUp;
  this.m_fov    = aFOV;
  this.m_near   = aNearClip;
  this.m_far    = aFarClip;
  this.m_width  = aSensorWidth;
  this.m_height = aSensorHeight;
}

/**
 * Get look at matrix 
 */
PerspectiveCamera.prototype.GetLookAtMatrix = function()
{
  return Matrix4.LookAt( this.m_pos , this.m_dir , this.m_up );
}

/*
* Get Projection matrix 
*/
PerspectiveCamera.prototype.GetProjectionMatrix = function()
{
  return Matrix4.Perspective( this.m_fov , this.m_width / this.m_height , this.m_near , this.m_far );
}

/**
 * Compute point a a plane passing through the destination point an parallel to the view plane 
 * The result is the intersection between the ray passing throug the eye to the pixel point and the plane  
 */
PerspectiveCamera.prototype.pointOnPlane = function( aPixX , aPixY )
{
  var diff = Vector.normalize( Vector.sub( this.m_dir , this.m_pos ) ); 

  // Plane passing through the destination point and facing the camera 
  var pl = Plane.initFromPositionAndNormal( this.m_dir , diff );
  // Point in the ray passing through the cam center and the pixel 
  var pt = this.UnProject( aPixX , aPixY );

  // Intersection between the plane and the line 
  return Plane.intersectLine( pl , this.m_pos , pt );
}

/**
 * Compute (nearest) point on a sphere centered at the destination point 
 * The result is the intersection between the ray passing throug the eye to the pixel point and the sphere 
 */
PerspectiveCamera.prototype.pointOnSphere = function( aPixX , aPixY , aRad )
{
  var pt = this.UnProject( aPixX , aPixY );
  var n = Vector.normalize( Vector.sub( this.m_pos , this.m_dir ) );

  // Plane passing through the destination point and facing the camera 
  var pl = Plane.initFromPositionAndNormal( this.m_dir , n );

  var hit_plane = Plane.intersectLine( pl , this.m_pos , pt );

  var z = hit_plane[2];

  var dst = Vector.sub( hit_plane , this.m_dir ); 

  var d = dst[0] * dst[0] + dst[1] * dst[1];

  if (d < aRad * Math.sqrt( 2 ) / 2.0 )
  {
    z = Math.sqrt( aRad * aRad - d * d );
  }
  else 
  {
    var t = aRad / Math.sqrt( 2 );
    z = t * t / d;
  }

  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = hit_plane[0] + n[0];
  res[1] = hit_plane[1] + n[1];
  res[2] = hit_plane[2] + n[2] * z; 

  return res; 
}

/**
 * Given a screen position, get a 3d point 
 */
PerspectiveCamera.prototype.UnProject = function( aPixX , aPixY )
{
  var view = this.GetLookAtMatrix( );
  var proj = this.GetProjectionMatrix( );

  var inv = Matrix4.inverse( Matrix4.transpose( Matrix4.mul( view , proj ) ) );

  var tmp = new Common.ARRAY_TYPE( 3 );
  tmp[ 0 ] = aPixX / this.m_width;
  tmp[ 1 ] = aPixY / this.m_height;
  tmp[ 2 ] = 0.0;

  var tmp2 = new Common.ARRAY_TYPE( 3 );
  tmp2[ 0 ] = 2.0 * tmp[ 0 ] - 1.0;
  tmp2[ 1 ] = 2.0 * tmp[ 1 ] - 1.0;
  tmp2[ 2 ] = 2.0 * tmp[ 2 ] - 1.0;

  return Matrix4.applyPoint( inv , tmp2 );
}

/**
 * Ensure sphere is inside view frustum 
 */
PerspectiveCamera.prototype.FitBoundingSphere = function( aCenter , aRadius )
{
  // Direction vector 
  var dir = Vector.normalize( Vector.sub( this.m_pos , this.m_dir ) ); 

  var arad = DegToRad( this.m_fov );
  var d = aRadius / Math.sin( arad * 0.5 );

  this.m_pos[0] = aCenter[0] + dir[0] * d; 
  this.m_pos[1] = aCenter[1] + dir[1] * d;
  this.m_pos[2] = aCenter[2] + dir[2] * d; 

  this.m_dir[0] = this.m_pos[0] - dir[0] * d; 
  this.m_dir[1] = this.m_pos[1] - dir[1] * d;
  this.m_dir[2] = this.m_pos[2] - dir[2] * d; 


}

PerspectiveCamera.prototype.zoom = function( aFactor )
{
  var dir = Vector.sub( this.m_dir , this.m_pos ); 
  var d = Vector.norm( dir );

  dir = Vector.normalize( dir ); 

  var new_pos = Vector.copy( this.m_pos ); 

  new_pos[0] += dir[0] * aFactor * d / 10.0; 
  new_pos[1] += dir[1] * aFactor * d / 10.0;
  new_pos[2] += dir[2] * aFactor * d / 10.0;

  var new_dir = Vector.normalize( Vector.sub( this.m_dir , new_pos ) );

  // Detect inversion of orientation 
  if (Vector.dot( new_dir , dir ) < 0.0 )
  {
    return; 
  }

  this.m_pos = new_pos; 
}

PerspectiveCamera.prototype.getPosition = function()
{
  return this.m_pos; 
}