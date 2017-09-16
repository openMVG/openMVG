// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @brief Constructor 
* @param x X-value 
* @param y Y-value 
* @param z Z-value 
* @param w W-value 
*/
Quaternion = function( aX , aY , aZ , aW )
{
  this.m_x = aX;
  this.m_y = aY;
  this.m_z = aZ;
  this.m_w = aW;
}

/**
 * Initialize a quaternion from an axis and an angle 
 */
Quaternion.prototype.setFromAxisAngle = function( anAxis , anAngle )
{
  var axis = Vector.normalize( anAxis );
  var sin_t = Math.sin( anAngle * 0.5 ); 

  if (Math.abs( anAngle ) > 0.0 )
  {
    this.m_w = Math.cos( anAngle * 0.5 );
    this.m_x = axis[0] * sin_t;
    this.m_y = axis[1] * sin_t;
    this.m_z = axis[2] * sin_t; 
  }
  else 
  {
    this.m_w = 0.0;
    this.m_x = 1.0;
    this.m_y = 0.0;
    this.m_z = 0.0; 
  }
}

/**
 * First imaginary part 
 */
Quaternion.prototype.getX = function( )
{
  return this.m_x;
}

/**
 * Second imaginary part 
 */
Quaternion.prototype.getY = function()
{
  return this.m_y;
}

/**
 * Third imaginary part 
 */
Quaternion.prototype.getZ = function()
{
  return this.m_z;
}

/**
 * Get Real part 
 */
Quaternion.prototype.getW = function()
{
  return this.m_w;
}

/**
 * Set first imaginary part 
 */
Quaternion.prototype.setX = function( aX )
{
  this.m_x = aX;
}

/**
 * Set second imaginary part 
 */
Quaternion.prototype.setY = function( aY )
{
  this.m_y = aY;
}

/**
 * Set third imaginary part 
 */
Quaternion.prototype.setZ = function( aZ )
{
  this.m_z = aZ;
}

/**
 * Set real part 
 */
Quaternion.prototype.setW = function( aW )
{
  this.m_w = aW;
}

/**
 * Conversion from 3x3 rotation matrix 
 * @note Rotation matrix is a 9 scalar vector 
 */
Quaternion.prototype.setFromMatrix3 = function( aMatrix )
{
  /**
   * 0 1 2 
   * 3 4 5
   * 6 7 8
   */
  this.m_w = Math.sqrt( 1.0 + aMatrix[ 0 ] + aMatrix[ 4 ] + aMatrix[ 8 ] ) / 2.0; 
  var fourW = 4.0 * this.m_w;
  this.m_x = ( aMatrix[ 7 ] - aMatrix[ 5 ] ) / fourW;
  this.m_y = ( aMatrix[ 2 ] - aMatrix[ 6 ] ) / fourW;
  this.m_z = ( aMatrix[ 3 ] - aMatrix[ 1 ] ) / fourW;
}

/**
 * Conversion from 3x3 rotation matrix 
 * @note Rotation matrix is a 9 scalar vector 
 */
Quaternion.prototype.setFromMatrix4 = function( aMatrix )
{
  /**
   * 0   1  2  3
   * 4   5  6  7
   * 8   9 10 11
   * 12 13 14 15
   */
  this.m_w = Math.sqrt( 1.0 + aMatrix[ 0 ] + aMatrix[ 5 ] + aMatrix[ 10 ] ) / 2.0; 
  var fourW = 4.0 * this.m_w;
  this.m_x = ( aMatrix[ 9 ] - aMatrix[ 6 ] ) / fourW;
  this.m_y = ( aMatrix[ 2 ] - aMatrix[ 8 ] ) / fourW;
  this.m_z = ( aMatrix[ 4 ] - aMatrix[ 1 ] ) / fourW;
}

/**
 * Get a 3x3 rotation matrix from this quaternion 
 */
Quaternion.prototype.toMatrix = function( )
{
  var res = Matrix3.create( );

  var sqw = this.m_w * this.m_w;
  var sqx = this.m_x * this.m_x;
  var sqy = this.m_y * this.m_y;
  var sqz = this.m_z * this.m_z; 

  var inv = 1.0 / ( sqx + sqy + sqz + sqw );

  // 0.0 -> 0
  res[0] = (   sqx - sqy - sqz + sqw ) * inv; 
  // 1.1 -> 0 
  res[4] = ( - sqx + sqy - sqz + sqw ) * inv;
  // 2.2 -> 0 
  res[8] = ( - sqx - sqy + sqz + sqw ) * inv; 

  var tmp1 = this.m_x * this.m_y;
  var tmp2 = this.m_z * this.m_w;

  // 1.0 -> 3
  res[3] = 2.0 * ( tmp1 + tmp2 ) * inv; 
  // 0.1 -> 1 
  res[1] = 2.0 * ( tmp1 - tmp2 ) * inv; 

  tmp1 = this.m_x * this.m_z;
  tmp2 = this.m_y * this.m_w; 

  // 2.0 -> 6 
  res[6] = 2.0 * ( tmp1 - tmp2 ) * inv; 
  // 0.2 -> 2 
  res[2] = 2.0 * ( tmp1 + tmp2 ) * inv; 

  tmp1 = this.m_y * this.m_z;
  tmp2 = this.m_x * this.m_w; 

  // 2.1 -> 7 
  res[7] = 2.0 * ( tmp1 + tmp2 ) * inv;
  // 1.2 -> 5 
  res[5] = 2.0 * ( tmp1 - tmp2 ) * inv; 

  return res; 
}


/**
 * Conjugate 
 */
Quaternion.prototype.conjugate = function( )
{
  return new Quaternion( -this.m_x , -this.m_y , -this.m_z , this.m_w );
}

/**
 * Norm 
 */
Quaternion.prototype.norm = function( )
{
  return Math.sqrt( this.m_x * this.m_x + 
                    this.m_y * this.m_y +
                    this.m_z * this.m_z +
                    this.m_w * this.m_w );
}

/**
 * Addition 
 */
Quaternion.prototype.add = function( anOther )
{
  return new Quaternion( this.m_x + anOther.m_x , 
                         this.m_y + anOther.m_y , 
                         this.m_z + anOther.m_z , 
                         this.m_w + anOther.m_w );
}

/**
 * Subtraction 
 */
Quaternion.prototype.sub = function( anOther )
{
  return new Quaternion( this.m_x - anOther.m_x , 
                         this.m_y - anOther.m_y ,
                         this.m_z - anOther.m_z ,
                         this.m_w - anOther.m_w ); 
}

/**
 * Product 
 */
Quaternion.prototype.mul = function( anOther )
{
  var w = this.m_w * anOther.m_w -
          this.m_x * anOther.m_x -
          this.m_y * anOther.m_y -
          this.m_z * anOther.m_z;

  var x = this.m_w * anOther.m_x +
          this.m_x * anOther.m_w +
          this.m_y * anOther.m_z -
          this.m_z * anOther.m_y;

  var y = this.m_w * anOther.m_y -
          this.m_x * anOther.m_z +
          this.m_y * anOther.m_w +
          this.m_z * anOther.m_x;

  var z = this.m_w * anOther.m_z +
          this.m_x * anOther.m_y -
          this.m_y * anOther.m_x + 
          this.m_z * anOther.m_w;  

  return new Quaternion( x , y , z , w );
}

/**
 * Scalar product 
 */
Quaternion.prototype.scalarProd = function( aVal )
{
  return new Quaternion( this.m_x * aVal , this.m_y * aVal , this.m_z * aVal , this.m_w * aVal );
}

/**
 * Inverse 
 */
Quaternion.prototype.inv = function( )
{
  var invnorm = 1.0 / ( this.m_x * this.m_x + 
                        this.m_y * this.m_y +
                        this.m_z * this.m_z +
                        this.m_w * this.m_w );
  return new Quaternion( -this.m_x * invnorm , 
                         -this.m_y * invnorm ,
                         -this.m_z * invnorm ,
                         this.m_w * invnorm );
}

/**
 * Compute unit quaternion 
 */
Quaternion.prototype.unit = function()
{
  var inv = 1.0 / Math.sqrt( this.m_x * this.m_x ,
                           this.m_y * this.m_y ,
                           this.m_z * this.m_z ,
                           this.m_w * this.m_w );
  return new Quaternion( this.m_x * inv ,
                         this.m_y * inv ,
                         this.m_z * inv ,
                         this.m_w * inv );
}
