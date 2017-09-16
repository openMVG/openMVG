// Copyright (c) 2016 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * Constructor 
 * @note identity translation (no rotation and translation)
 */
DualQuaternion = function()
{
  this.m_r = new Quaternion( 0.0 , 0.0 , 0.0 , 1.0 );
  this.m_d = new Quaternion( 0.0 , 0.0 , 0.0 , 0.0 );
}

/**
 * Initialize a dual quaternion from a quaternion rotation 
 */
DualQuaternion.prototype.setFromRotationQuaternion = function( aQuat )
{
  this.m_r = new Quaternion( aQuat.getX() , aQuat.getY() , aQuat.getZ() , aQuat.getW() );
  this.m_d = new Quaternion( 0.0 , 0.0 , 0.0 , 0.0 );
}

/**
 * Initialize a dual quaternion from a translation vector 
 */
DualQuaternion.prototype.setFromTranslationVector = function( aVec )
{
  this.m_r = new Quaternion( 0.0 , 0.0 , 0.0 , 1.0 );
  this.m_d = new Quaternion( aVec[0] * 0.5 , aVec[1] * 0.5 , aVec[2] * 0.5 , 0.0 );
}

/**
 * Deep copy of a dual quaternion 
 */
DualQuaternion.prototype.copy = function( aDualQuat )
{
  this.m_r = new Quaternion( aDualQuat.m_r.getX() , aDualQuat.m_r.getY() , aDualQuat.m_r.getZ() , aDualQuat.m_r.getW() );
  this.m_d = new Quaternion( aDualQuat.m_d.getX() , aDualQuat.m_d.getY() , aDualQuat.m_d.getZ() , aDualQuat.m_d.getW() );
}

/**
 * Product of two dual quaternion 
 */
DualQuaternion.prototype.mul= function( aDualQuat )
{
  var res = new DualQuaternion();

  res.m_r = this.m_r.mul( aDualQuat.m_r );
  var tmp1 = this.m_r.mul( aDualQuat.m_d );
  var tmp2 = this.m_d.mul( aDualQuat.m_r );
  res.m_d = tmp1.add( tmp2 ); 

  return res; 
}

/**
 * Convert dual quaternion to a 4x4 matrix 
 */
DualQuaternion.prototype.toMatrix = function( )
{
  var res = new Array( 16 ); 

  // Rotational part (a 3x3 matrix)
  var rot = this.m_r.toMatrix(); 

  /**
   * 0   1   2   3 
   * 4   5   6   7
   * 8   9  10  11
   * 12 13  14  15
   */

  res[0] = rot[0];
  res[1] = rot[1];
  res[2] = rot[2];
  
  res[4] = rot[3];
  res[5] = rot[4];
  res[6] = rot[5];

  res[8] = rot[6];
  res[9] = rot[7];
  res[10] = rot[8]; 

  // Translational part 
  var q = ( this.m_d.scalarProd( 2.0 ) ).mul( this.m_r.conjugate() );

  res[3]  = q.getX();
  res[7]  = q.getY();
  res[11] = q.getZ(); 

  res[12] = 0.0;
  res[13] = 0.0;
  res[14] = 0.0;
  res[15] = 1.0;  

  return res; 
}
