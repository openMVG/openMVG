// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Class creating a plane */
var Plane = {};

/* create a plane using it's normal and position */
Plane.initFromPositionAndNormal = function( aPos , aNormal )
{
  var res = new Common.ARRAY_TYPE( 6 );

  res[ 0 ] = aNormal[ 0 ];
  res[ 1 ] = aNormal[ 1 ];
  res[ 2 ] = aNormal[ 2 ];
  res[ 3 ] = aPos[ 0 ];
  res[ 4 ] = aPos[ 1 ];
  res[ 5 ] = aPos[ 2 ];

  return res;
}

/* point of intersection between a plane and a line passing through two points */
Plane.intersectLine = function( aPlane , aPt1 , aPt2 )
{
  var l = Vector.sub( aPt2 , aPt1 ); 
  var n = Vector.create( aPlane[0] , aPlane[1] , aPlane[2] ); 
  var p = Vector.create( aPlane[3] , aPlane[4] , aPlane[5] ); 
  var diff = Vector.sub( p , aPt1 ); 

  var d = Vector.dot( diff , n ) / Vector.dot( l , n );

  var res = new Common.ARRAY_TYPE( 3 );

  res[ 0 ] = aPt1[ 0 ] + d * l[ 0 ];
  res[ 1 ] = aPt1[ 1 ] + d * l[ 1 ]; 
  res[ 2 ] = aPt1[ 2 ] + d * l[ 2 ];

  return res;
}
