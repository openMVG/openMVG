// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * Global object used to manage a 3x3 matrix 
 */
var Matrix3 = {} 

/**
 * Create an identity matrix 
 */
Matrix3.create = function()
{
  var res = new Common.ARRAY_TYPE( 9 );

  res[0] = 1.0;
  res[1] = 0.0;
  res[2] = 0.0;

  res[3] = 0.0;
  res[4] = 1.0;
  res[5] = 0.0;

  res[6] = 0.0;
  res[7] = 0.0;
  res[8] = 1.0;

  return res; 
}

/**
 * Deep copy 
 */
Matrix3.copy = function( aMatrix )
{
  var res = new Common.ARRAY_TYPE( 9 );

  res[0] = aMatrix[0];
  res[1] = aMatrix[1];
  res[2] = aMatrix[2];

  res[3] = aMatrix[3];
  res[4] = aMatrix[4];
  res[5] = aMatrix[5];

  res[6] = aMatrix[6];
  res[7] = aMatrix[7];
  res[8] = aMatrix[8];

  return out; 
}

/**
 * Addition of two matrices 
 */
Matrix3.add = function( aMatrix1 , aMatrix2 )
{
  var res = new Common.ARRAY_TYPE( 9 );

  res[0] = aMatrix1[0] + aMatrix2[0];
  res[1] = aMatrix1[1] + aMatrix2[1];
  res[2] = aMatrix1[2] + aMatrix2[2];

  res[3] = aMatrix1[3] + aMatrix2[3];
  res[4] = aMatrix1[4] + aMatrix2[4];
  res[5] = aMatrix1[5] + aMatrix2[5];

  res[6] = aMatrix1[6] + aMatrix2[6];
  res[7] = aMatrix1[7] + aMatrix2[7];
  res[8] = aMatrix1[8] + aMatrix2[8];

  return res; 
}

/**
 * Subtraction of two matrices 
 */
Matrix3.sub = function( aMatrix1 , aMatrix2 ) 
{
  var res = new Common.ARRAY_TYPE( 9 );

  res[0] = aMatrix1[0] - aMatrix2[0];
  res[1] = aMatrix1[1] - aMatrix2[1];
  res[2] = aMatrix1[2] - aMatrix2[2];

  res[3] = aMatrix1[3] - aMatrix2[3];
  res[4] = aMatrix1[4] - aMatrix2[4];
  res[5] = aMatrix1[5] - aMatrix2[5];

  res[6] = aMatrix1[6] - aMatrix2[6];
  res[7] = aMatrix1[7] - aMatrix2[7];
  res[8] = aMatrix1[8] - aMatrix2[8];

  return res; 
}

/**
 * Product of two matrices 
 */
Matrix3.mul = function( aMatrix1 , aMatrix2 )
{
  var res = new Common.ARRAY_TYPE( 9 );

  /**
   * Standard product : 
   *                 b0 b1 b2
   *                 b3 b4 b5
   *                 b6 b7 b8
   *                
   *      a0 a1 a2   c0 c1 c2
   *      a3 a4 a5   c3 c4 c5
   *      a6 a7 a8   c6 c7 c8
   * 
   * 
   * c0 =  a0 * b0 +  a1 * b3 +  a2 * b6
   * c1 =  a0 * b1 +  a1 * b4 +  a2 * b7
   * ... 
   */

  res[0] = aMatrix1[0] * aMatrix2[0] + aMatrix1[1] * aMatrix2[3] + aMatrix1[2] * aMatrix2[6];
  res[1] = aMatrix1[0] * aMatrix2[1] + aMatrix1[1] * aMatrix2[4] + aMatrix1[2] * aMatrix2[7];
  res[2] = aMatrix1[0] * aMatrix2[2] + aMatrix1[1] * aMatrix2[5] + aMatrix1[2] * aMatrix2[8];

  res[3] = aMatrix1[3] * aMatrix2[0] + aMatrix1[4] * aMatrix2[3] + aMatrix1[5] * aMatrix2[6];
  res[4] = aMatrix1[3] * aMatrix2[1] + aMatrix1[4] * aMatrix2[4] + aMatrix1[5] * aMatrix2[7];
  res[5] = aMatrix1[3] * aMatrix2[2] + aMatrix1[4] * aMatrix2[5] + aMatrix1[5] * aMatrix2[8];

  res[6] = aMatrix1[6] * aMatrix2[0] + aMatrix1[7] * aMatrix2[3] + aMatrix1[8] * aMatrix2[6];
  res[7] = aMatrix1[6] * aMatrix2[1] + aMatrix1[7] * aMatrix2[4] + aMatrix1[8] * aMatrix2[7];
  res[8] = aMatrix1[6] * aMatrix2[2] + aMatrix1[7] * aMatrix2[5] + aMatrix1[8] * aMatrix2[8];

  return res; 
}

/**
 * Transposition of a matrix 
 */
Matrix3.transpose = function( aMatrix )
{
  var res = new Common.ARRAY_TYPE( 9 );

  res[0] = aMatrix[0];
  res[1] = aMatrix[3];
  res[2] = aMatrix[6];
  
  res[3] = aMatrix[1];
  res[4] = aMatrix[4];
  res[5] = aMatrix[7];

  res[6] = aMatrix[2];
  res[7] = aMatrix[5];
  res[8] = aMatrix[8];

  return res; 
}
