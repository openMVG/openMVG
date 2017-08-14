// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

var Matrix4 = {};

/**
 * Create an identity matrix 
 */
Matrix4.create = function()
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = 1.0;
  res[1] = 0.0;
  res[2] = 0.0;
  res[3] = 0.0;

  res[4] = 0.0;
  res[5] = 1.0;
  res[6] = 0.0;
  res[7] = 0.0;

  res[8]  = 0.0;
  res[9]  = 0.0;
  res[10] = 1.0; 
  res[11] = 0.0;

  res[12] = 0.0;
  res[13] = 0.0;
  res[14] = 0.0;
  res[15] = 1.0;

  return res; 
}

/* Copy a matrix */
Matrix4.copy = function( aMatrix )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = aMatrix[0];
  res[1] = aMatrix[1];
  res[2] = aMatrix[2];
  res[3] = aMatrix[3];

  res[4] = aMatrix[4];
  res[5] = aMatrix[5];
  res[6] = aMatrix[6];
  res[7] = aMatrix[7];

  res[8] = aMatrix[8];
  res[9] = aMatrix[9];
  res[10] = aMatrix[10];
  res[11] = aMatrix[11];

  res[12] = aMatrix[12];
  res[13] = aMatrix[13];
  res[14] = aMatrix[14];
  res[15] = aMatrix[15];

  return res;   
}

/**
 * Sum of two matrices 
 */
Matrix4.add = function( aMatrix1 , aMatrix2 )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = aMatrix1[0] + aMatrix2[0];
  res[1] = aMatrix1[1] + aMatrix2[1];
  res[2] = aMatrix1[2] + aMatrix2[2];
  res[3] = aMatrix1[3] + aMatrix2[3];

  res[4] = aMatrix1[4] + aMatrix2[4];
  res[5] = aMatrix1[5] + aMatrix2[5];
  res[6] = aMatrix1[6] + aMatrix2[6];
  res[7] = aMatrix1[7] + aMatrix2[7];

  res[8] = aMatrix1[8] + aMatrix2[8];
  res[9] = aMatrix1[9] + aMatrix2[9];
  res[10] = aMatrix1[10] + aMatrix2[10];
  res[11] = aMatrix1[11] + aMatrix2[11];

  res[12] = aMatrix1[12] + aMatrix2[12];
  res[13] = aMatrix1[13] + aMatrix2[13];
  res[14] = aMatrix1[14] + aMatrix2[14];
  res[15] = aMatrix1[15] + aMatrix2[15];

  return res;   
}

/**
 * Subtraction of two matrices 
 */
Matrix4.sub = function( aMatrix1 , aMatrix2 )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = aMatrix1[0] - aMatrix2[0];
  res[1] = aMatrix1[1] - aMatrix2[1];
  res[2] = aMatrix1[2] - aMatrix2[2];
  res[3] = aMatrix1[3] - aMatrix2[3];

  res[4] = aMatrix1[4] - aMatrix2[4];
  res[5] = aMatrix1[5] - aMatrix2[5];
  res[6] = aMatrix1[6] - aMatrix2[6];
  res[7] = aMatrix1[7] - aMatrix2[7];

  res[8] = aMatrix1[8] - aMatrix2[8];
  res[9] = aMatrix1[9] - aMatrix2[9];
  res[10] = aMatrix1[10] - aMatrix2[10];
  res[11] = aMatrix1[11] - aMatrix2[11];

  res[12] = aMatrix1[12] - aMatrix2[12];
  res[13] = aMatrix1[13] - aMatrix2[13];
  res[14] = aMatrix1[14] - aMatrix2[14];
  res[15] = aMatrix1[15] - aMatrix2[15];

  return res;     
}

/**
 * Matrix product 
 */
Matrix4.mul = function( a , b )
{
  var res = new Common.ARRAY_TYPE( 16 );

  /**
   * Standard matrix multiplication 
   * 
   * 
   *                     b0  b1  b2  b3 
   *                     b4  b5  b6  b7
   *                     b8  b9 b10 b11
   *                    b12 b13 b14 b15 
   * 
   *   a0  a1  a2  a3    c0  c1  c2  c3
   *   a4  a5  a6  a7    c4  c5  c6  c7
   *   a8  a9 a10 a11    c8  c9 c10 c11
   *  a12 a13 a14 a15   c12 c13 c14 c15
   * 
   *  c0 = a0 * b0  +  a1 * b4  +  a2 * b8  +  a3 * b12
   *  c1 = ... 
   */

  res[0]  = a[0] * b[0]  + a[1] * b[4]   + a[2] * b[8]    + a[3] * b[12];
  res[1]  = a[0] * b[1]  + a[1] * b[5]   + a[2] * b[9]    + a[3] * b[13];
  res[2]  = a[0] * b[2]  + a[1] * b[6]   + a[2] * b[10]   + a[3] * b[14];
  res[3]  = a[0] * b[3]  + a[1] * b[7]   + a[2] * b[11]   + a[3] * b[15];
  
  res[4]  = a[4] * b[0]  + a[5] * b[4]   + a[6] * b[8]    + a[7] * b[12];
  res[5]  = a[4] * b[1]  + a[5] * b[5]   + a[6] * b[9]    + a[7] * b[13];
  res[6]  = a[4] * b[2]  + a[5] * b[6]   + a[6] * b[10]   + a[7] * b[14];
  res[7]  = a[4] * b[3]  + a[5] * b[7]   + a[6] * b[11]   + a[7] * b[15];

  res[8]  = a[8] * b[0]  + a[9] * b[4]   + a[10] * b[8]   + a[11] * b[12];
  res[9]  = a[8] * b[1]  + a[9] * b[5]   + a[10] * b[9]   + a[11] * b[13];
  res[10] = a[8] * b[2]  + a[9] * b[6]   + a[10] * b[10]  + a[11] * b[14];
  res[11] = a[8] * b[3]  + a[9] * b[7]   + a[10] * b[11]  + a[11] * b[15];

  res[12] = a[12] * b[0] + a[13] * b[4]  + a[14] * b[8]   + a[15] * b[12];
  res[13] = a[12] * b[1] + a[13] * b[5]  + a[14] * b[9]   + a[15] * b[13];
  res[14] = a[12] * b[2] + a[13] * b[6]  + a[14] * b[10]  + a[15] * b[14];
  res[15] = a[12] * b[3] + a[13] * b[7]  + a[14] * b[11]  + a[15] * b[15];

  return res; 
}

/**
 * Transposition of a matrix 
 */
Matrix4.transpose = function( aMatrix )
{
  var res = new Common.ARRAY_TYPE( 16 );

  /**
   * 0   1  2  3 
   * 4   5  6  7
   * 8   9 10 11
   * 12 13 14 15
   */

  res[0] = aMatrix[0];
  res[1] = aMatrix[4];
  res[2] = aMatrix[8];
  res[3] = aMatrix[12];

  res[4] = aMatrix[1];
  res[5] = aMatrix[5];
  res[6] = aMatrix[9];
  res[7] = aMatrix[13];

  res[8] = aMatrix[2];
  res[9] = aMatrix[6];
  res[10] = aMatrix[10];
  res[11] = aMatrix[14];

  res[12] = aMatrix[3];
  res[13] = aMatrix[7];
  res[14] = aMatrix[11];
  res[15] = aMatrix[15];

  return res; 
}

/**
 * Inverse of a matrix 
 */
Matrix4.inverse = function( aMatrix )
{
  /* see intel paper: streaming SIMD Extensions - Inverse of 4x4 matrix   */
  var  tmp = new Common.ARRAY_TYPE( 16 ); /* temp array for pairs             */
  var  src = new Common.ARRAY_TYPE( 16 ); /* array of transpose source matrix */
  var out = new Common.ARRAY_TYPE( 16 );

  /* transpose matrix */
  for (var i = 0; i < 4; ++i)
  {
    src[i]        = aMatrix[i*4];
    src[i + 4]    = aMatrix[i*4 + 1];
    src[i + 8]    = aMatrix[i*4 + 2];
    src[i + 12]   = aMatrix[i*4 + 3];
  }

/* calculate pairs for first 8 elements (cofactors) */
  tmp[0]  = src[10] * src[15];
  tmp[1]  = src[11] * src[14];
  tmp[2]  = src[9]  * src[15];
  tmp[3]  = src[11] * src[13];
  tmp[4]  = src[9]  * src[14];
  tmp[5]  = src[10] * src[13];
  tmp[6]  = src[8]  * src[15];
  tmp[7]  = src[11] * src[12];
  tmp[8]  = src[8]  * src[14];
  tmp[9]  = src[10] * src[12];
  tmp[10] = src[8]  * src[13];
  tmp[11] = src[9]  * src[12];
  /* calculate first 8 elements (cofactors) */
  out[0]  = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
  out[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  out[1]  = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
  out[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
  out[2]  = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  out[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
  out[3]  = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
  out[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
  out[4]  = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
  out[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
  out[5]  = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
  out[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  out[6]  = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
  out[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  out[7]  = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  out[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0]  = src[2]*src[7];
  tmp[1]  = src[3]*src[6];
  tmp[2]  = src[1]*src[7];
  tmp[3]  = src[3]*src[5];
  tmp[4]  = src[1]*src[6];
  tmp[5]  = src[2]*src[5];
  tmp[6]  = src[0]*src[7];
  tmp[7]  = src[3]*src[4];
  tmp[8]  = src[0]*src[6];
  tmp[9]  = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  out[8]  = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
  out[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  out[9]  = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
  out[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
  out[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  out[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
  out[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  out[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
  out[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  out[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
  out[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
  out[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
  out[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
  out[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
  out[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
  out[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  var det=src[0]*out[0]+src[1]*out[1]+src[2]*out[2]+src[3]*out[3];
  /* calculate matrix inverse */
  det = 1.0 / det;
  for (var j = 0; j < 16; ++j )
   out[j] *= det;

  return out;
}

/**
 * Create a translation matrix 
 */
Matrix4.createTranslation = function( aVec )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = 1.0;
  res[1] = 0.0; 
  res[2] = 0.0;
  res[3] = aVec[0];

  res[4] = 0.0;
  res[5] = 1.0; 
  res[6] = 0.0;
  res[7] = aVec[1];

  res[8] = 0.0;
  res[9] = 0.0;
  res[10] = 1.0;
  res[11] = aVec[2];

  res[12] = 0.0;
  res[13] = 0.0;
  res[14] = 0.0;
  res[15] = 1.0;

  return res; 
}

/**
 * Create a scale matrix 
 */
Matrix4.createScale = function( aScaleX , aScaleY , aScaleZ )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = aScaleX;
  res[1] = 0.0; 
  res[2] = 0.0;
  res[3] = 0.0;

  res[4] = 0.0;
  res[5] = aScaleY;
  res[6] = 0.0;
  res[7] = 0.0; 

  res[8] = 0.0;
  res[9] = 0.0;
  res[10] = aScaleZ;
  res[11] = 0.0;

  res[12] = 0.0;
  res[13] = 0.0;
  res[14] = 0.0;
  res[15] = 1.0; 

  return res; 
}

/**
 * Create a uniform scaling matrix 
 */
Matrix4.createUniformScale = function( aScale )
{
  var res = new Common.ARRAY_TYPE( 16 );

  res[0] = aScale;
  res[1] = 0.0; 
  res[2] = 0.0;
  res[3] = 0.0;

  res[4] = 0.0;
  res[5] = aScale;
  res[6] = 0.0;
  res[7] = 0.0; 

  res[8] = 0.0;
  res[9] = 0.0;
  res[10] = aScale;
  res[11] = 0.0;

  res[12] = 0.0;
  res[13] = 0.0;
  res[14] = 0.0;
  res[15] = 1.0; 

  return res; 
}

/**
 * Apply current geometric transformation to a point 
 */
Matrix4.applyPoint = function( aMat , aPoint )
{
  // Matrix product 
  var x = aMat[0] * aPoint[0] + aMat[1] * aPoint[1] + aMat[2] * aPoint[2] + aMat[3];
  var y = aMat[4] * aPoint[0] + aMat[5] * aPoint[1] + aMat[6] * aPoint[2] + aMat[7];
  var z = aMat[8] * aPoint[0] + aMat[9] * aPoint[1] + aMat[10] * aPoint[2] + aMat[11];

  var w = 1.0 / ( aMat[12] * aPoint[0] + aMat[13] * aPoint[1] + aMat[14] * aPoint[2] + aMat[15] );

  var res = new Common.ARRAY_TYPE( 3 );

  res[0] = x * w;
  res[1] = y * w;
  res[2] = z * w; 

  return res; 
}

/**
 * Look at matrix 
 */
Matrix4.LookAt = function( aFrom , aDest , aUp )
{
  var eye = aFrom;
  var dir = Vector.normalize( Vector.sub( aFrom , aDest ) ); 
  var up = Vector.normalize( aUp );
  var right = Vector.normalize( Vector.cross( up , dir ) );
  up = Vector.cross( dir , right );


  var m = new Common.ARRAY_TYPE( 16 ); 

      m[0] = right[0];
      m[1] = up[0];
      m[2] = dir[0];
      m[3] = 0.0;

      m[4] = right[1];
      m[5] = up[1];
      m[6] = dir[1];
      m[7] = 0.0;

      m[8] = right[2];
      m[9] = up[2];
      m[10] = dir[2];
      m[11] = 0.0;

      m[12] = - Vector.dot( right , eye );
      m[13] = - Vector.dot( up , eye );
      m[14] = - Vector.dot( dir , eye );
      m[15] = 1.0;

      return m; 
}

/**
 * Perspective matrix 
*/
Matrix4.Perspective = function( aFov , aPixelRatio , aNear , aFar ) 
{

      var range  = Math.tan( aFov * Math.PI / 360.0 ) * aNear;
      var left   = -range * aPixelRatio;
      var right  = range * aPixelRatio;
      var bottom = - range;
      var top    = range;

      var m = new Common.ARRAY_TYPE( 16 ); 

      m[0] =  ( 2.0 * aNear ) / ( right - left );
      m[1] = 0.0;
      m[2] = 0.0;
      m[3] = 0.0;

      m[4] = 0.0;
      m[5] = ( 2.0 * aNear ) / ( top - bottom );
      m[6] = 0.0;
      m[7] = 0.0;

      m[8] = ( right + left ) / ( right - left );
      m[9] = ( top + bottom ) / ( top - bottom );
      m[10] = - ( aFar + aNear ) / ( aFar - aNear );
      m[11] = -1.0;

      m[12] = 0.0;
      m[13] = 0.0;
      m[14] = - ( 2.0 * aFar * aNear ) / ( aFar - aNear );
      m[15] = 0.0;

      return m;
}