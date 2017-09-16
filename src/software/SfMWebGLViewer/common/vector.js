// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* A simple 3d vector class */
Vector = {} 

/* Create a vector by defining all it's components */
Vector.create = function( aX , aY , aZ ) 
{
  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = aX;
  res[1] = aY;
  res[2] = aZ; 

  return res; 
}

/* Copy a vector */
Vector.copy = function( aVec )
{
  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = aVec[0];
  res[1] = aVec[1];
  res[2] = aVec[2];

  return res; 
}

/* Add two vectors */
Vector.add = function( aVector1 , aVector2 )
{
  var res = new Common.ARRAY_TYPE( 3 ); 

  res[0] = aVector1[0] + aVector2[0];
  res[1] = aVector1[1] + aVector2[1];
  res[2] = aVector1[2] + aVector2[2];

  return res;   
}

/* Subtract two vectors */
Vector.sub = function( aVector1 , aVector2 )
{
  var res = new Common.ARRAY_TYPE( 3 ); 

  res[0] = aVector1[0] - aVector2[0];
  res[1] = aVector1[1] - aVector2[1];
  res[2] = aVector1[2] - aVector2[2];

  return res; 
}

/* Dot product of two vectors */
Vector.dot = function( aVector1 , aVector2 ) 
{
  return aVector1[0] * aVector2[0] + 
         aVector1[1] * aVector2[1] +
         aVector1[2] * aVector2[2];
}

/* Norm of a vector */
Vector.norm = function( aVector )
{
  return Math.sqrt( aVector[0] * aVector[0] + aVector[1] * aVector[1] + aVector[2] * aVector[2] );
}

/* Normalize a vector */
Vector.normalize = function( aVector )
{
  var inv = 1.0 / Math.sqrt( aVector[0] * aVector[0] + aVector[1] * aVector[1] + aVector[2] * aVector[2] );

  return Vector.create( aVector[0] * inv , aVector[1] * inv , aVector[2] * inv ); 
}

/* Cross product of two vectors */
Vector.cross = function( aVector1 , aVector2 )
{
  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = aVector1[1] * aVector2[2] - aVector1[2] * aVector2[1];
  res[1] = aVector1[2] * aVector2[0] - aVector1[0] * aVector2[2];
  res[2] = aVector1[0] * aVector2[1] - aVector1[1] * aVector2[0];

  return res;  
}

/* Negation of a vector */
Vector.negate = function( aVector )
{
  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = - aVector[0];
  res[1] = - aVector[1];
  res[2] = - aVector[2];

  return res; 
}

