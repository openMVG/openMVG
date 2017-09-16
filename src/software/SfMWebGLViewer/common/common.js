// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
* Define a global object containing various common definitions 
*/
var Common = {} 

// The standard type used to create a floating point array 
Common.ARRAY_TYPE =  (typeof Float32Array !== 'undefined') ? Float32Array : Array;
Common.INT_ARRAY_TYPE = (typeof Int32Array !== 'undefined' ) ? Int32Array : Array;

/** 
 * Convert an angle from degree to radian 
 */
function DegToRad( aDeg ) 
{
  return aDeg * 0.01745329252;
}

/**
 * Convert an angle from radian to degree 
 */
function RadToDeg( aRad )
{
  return aRad * 57.295779513;
}


/* Check if DOM element has a class property */
function HasClass( aNode , aClassName )
{
  if (aNode.className )
  {
    return aNode.className.match( new RegExp( '(\\s|^)' + aClassName + '(\\s|$)') );
  }
  else
  {
    return false;
  }
}

/* Add a classname to a DOM element */
function AddClass( aNode , aClassName )
{
  aNode.className += " " + aClassName;
}

/* Remove a classname to a DOM element */
function RemoveClass( aNode , aClassName )
{
  if (HasClass( aNode , aClassName ) )
  {
    var reg = new RegExp('(\\s|^)' + aClassName + '(\\s|$)');
    aNode.className = aNode.className.replace(reg, ' ');
  }
}

/* Compute bounding sphere of a point list 
* Return a 4d array (cx,cy,cz,rad) */
function ComputeBoundingSphere( aPointList )
{
  var min_x = 10e10;
  var min_y = 10e10;
  var min_z = 10e10;
  var max_x = -10e10;
  var max_y = -10e10;
  var max_z = -10e10;

  for (var i = 0; i < aPointList.length / 3; ++i )
  {
    min_x = aPointList[ 3 * i ] < min_x     ? aPointList[ 3 * i ]     : min_x;
    min_y = aPointList[ 3 * i + 1 ] < min_y ? aPointList[ 3 * i + 1 ] : min_y;
    min_z = aPointList[ 3 * i + 2 ] < min_z ? aPointList[ 3 * i + 2 ] : min_z;

    max_x = aPointList[ 3 * i ] > max_x ? aPointList[ 3 * i ] : max_x;
    max_y = aPointList[ 3 * i + 1 ] > max_y ? aPointList[ 3 * i + 1 ] : max_y;
    max_z = aPointList[ 3 * i + 2 ] > max_z ? aPointList[ 3 * i + 2 ] : max_z;    
  }

  var res = new Common.ARRAY_TYPE( 4 );
  res[0] = ( min_x + max_x ) / 2.0;
  res[1] = ( min_y + max_y ) / 2.0;
  res[2] = ( min_z + max_z ) / 2.0; 
  var dx = max_x - res[0];
  var dy = max_y - res[1];
  var dz = max_z - res[2]; 
  res[3] = Math.sqrt( dx * dx + dy * dy + dz * dz );
  return res; 
}