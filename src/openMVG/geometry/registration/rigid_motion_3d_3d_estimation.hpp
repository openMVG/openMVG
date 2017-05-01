// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_REGISTRATION_RIGID_MOTION_3D_3D_ESTIMATION_HPP
#define OPENMVG_GEOMETRY_REGISTRATION_RIGID_MOTION_3D_3D_ESTIMATION_HPP

#include "openMVG/numeric/numeric.h"

#include <ceres/ceres.h>

namespace openMVG
{
namespace geometry
{
namespace registration
{

/**
* @brief Function used to compute error using point to point distance metric
*/
struct RigidMotion3d3dEstimationCost
{
  public:

    /**
    * @brief Constructor
    * @param ref Reference point
    * @param data Data point (ie : point to transform)
    */
    RigidMotion3d3dEstimationCost( const Vec3 &ref, const Vec3 &data )
      : m_ref( ref ),
        m_data( data )
    {
    }

    /**
    * @brief Compute distance error after applying a transformation to a point
    * @param quat A quaternion
    * @param tra A translation
    * @param[out] residual
    * @note quat is a 4 component array
    * @note tra is a 3 component array
    * @note residual is a 3 component array (residual on x,y,z)
    */
    template<typename T>
    bool operator()( const T * const quat , const T * const tra , T * residual ) const
    {
      // The points (promote to ceres::jet type )
      Eigen::Matrix<T, 3, 1> point;
      point << T( m_data[0] ), T( m_data[1] ), T( m_data[2] );

      Eigen::Matrix<T, 3, 1> ref;
      ref << T( m_ref[ 0 ] ), T( m_ref[ 1 ] ), T( m_ref[ 2 ] );

      // Map params to Eigen transformations
      const Eigen::Quaternion<T> q = Eigen::Map<const Eigen::Quaternion<T> >( quat );
      const Eigen::Matrix<T, 3, 1> t = Eigen::Map<const Eigen::Matrix<T, 3, 1> >( tra );

      // Compute transformation
      const Eigen::Matrix<T, 3, 1> p = ( q * point ) + t ;

      residual[ 0 ] = p[0] - ref[0] ;
      residual[ 1 ] = p[1] - ref[1] ;
      residual[ 2 ] = p[2] - ref[2] ;

      return true ;
    }

  private:
    /// The reference point
    Vec3 m_ref;
    /// The data point
    Vec3 m_data;
};

/**
* @brief Function used to compute error using point to plane distance metric
*/
struct RigidMotion3d3dNormalEstimationCost
{
  public:

    /**
    * @brief Constructor
    * @param ref Reference point
    * @param ref_n Normal at reference point
    * @param data Data point (ie : point to transform)
    */
    RigidMotion3d3dNormalEstimationCost( const Vec3 &ref, const Vec3 & ref_n , const Vec3 &data )
      : m_ref( ref ),
        m_ref_n( ref_n ) ,
        m_data( data )
    {
    }


    /**
    * @brief Compute distance error after applying a transformation to a point
    * @param quat A quaternion
    * @param tra A translation
    * @param[out] residual
    * @note quat is a 4 component array
    * @note tra is a 3 component array
    * @note residual is a 1 component array (residual on distance to plane)
    */
    template<typename T>
    bool operator()( const T * const quat , const T * const tra , T * residual ) const
    {
      // The points
      Eigen::Matrix<T, 3, 1> point;
      point << T( m_data[0] ), T( m_data[1] ), T( m_data[2] );

      Eigen::Matrix<T, 3, 1> ref;
      ref << T( m_ref[ 0 ] ), T( m_ref[ 1 ] ), T( m_ref[ 2 ] );

      Eigen::Matrix<T, 3, 1> ref_n;
      ref_n << T( m_ref_n[0] ) , T( m_ref_n[1] ) , T( m_ref_n[2] ) ;

      // The transformation
      const Eigen::Quaternion<T> q = Eigen::Map<const Eigen::Quaternion<T> >( quat );
      const Eigen::Matrix<T, 3, 1> t = Eigen::Map<const Eigen::Matrix<T, 3, 1> >( tra );

      // Compute transformation
      const Eigen::Matrix<T, 3, 1> p = ( q * point ) + t ;

      residual[ 0 ] = ( p - ref ).dot( ref_n ) ;

      return true ;
    }

  private:
    Vec3 m_ref ;
    Vec3 m_ref_n ;
    Vec3 m_data ;
};

/**
* @brief Ceres interface used to compute rigid transformation between two set of 3d points
* @note this use either the point to point distance metric or point to plane
*/
template <typename Scalar>
class RigidMotion3d3dEstimation
{
    using point_list_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>;

  public:
    RigidMotion3d3dEstimation();

    /**
    * @brief Computes ridig motion (q,t) : q * data + t that maps data on ref
    * @param ref Reference point cloud
    * @param data Data point cloud
    * @param corresps (id of the correspondances from data to ref, -1 if no corresp)
    * @return R,t Rigid motion that minimise the MSE
    */
    std::pair<Eigen::Quaternion<Scalar>, Eigen::Matrix<Scalar, 3, 1>> operator()( const point_list_type &ref, const point_list_type &data, const std::vector<int> &corresps ) const ;

    /**
    * @brief Computes ridig motion (q,t) : q * data + t that maps data on ref
    * @param ref Reference point cloud
    * @param ref_normal Reference point cloud normals
    * @param data Data point cloud
    * @param corresps (id of the correspondances from data to ref, -1 if no corresp)
    * @return < quat , tra > Rigid motion that minimise the MSE of the point to plane distance
    */
    std::pair<Eigen::Quaternion<Scalar>, Eigen::Matrix<Scalar, 3, 1>> operator()( const point_list_type &ref,
        const point_list_type &ref_normal,
        const point_list_type &data,
        const std::vector<int> &corresps ) const ;


  private:
};

template <typename Scalar>
RigidMotion3d3dEstimation<Scalar>::RigidMotion3d3dEstimation()
{
}

/**
  * @brief Computes ridig motion (R,t) : R * data + t that maps data on ref
  * @param ref Reference point cloud
  * @param data Data point cloud
  * @param corresps (id of the correspondances from data to ref, -1 if no corresp)
  * @return R,t Rigid motion that minimise the MSE
  */
template <typename Scalar>
std::pair<Eigen::Quaternion<Scalar>, Eigen::Matrix<Scalar, 3, 1>> RigidMotion3d3dEstimation<Scalar>::operator()( const point_list_type &ref,
    const point_list_type &data,
    const std::vector<int> &corresps ) const
{
  // 1 Add residuals functions
  ceres::Problem problem;
  // Initial value (no transformation)
  Eigen::Quaternion<Scalar> q = Eigen::Quaternion<Scalar>::Identity();
  Eigen::Matrix<Scalar, 3, 1> t( 0, 0, 0 );

  for ( size_t id_data = 0; id_data < data.rows() ; ++id_data )
  {
    if ( corresps[ id_data ] >= 0 )
    {
      const int id_ref  = corresps[ id_data ];

      const Vec3 ref_pt = ref.row( id_ref ).template cast<double>() ;
      const Vec3 data_pt = data.row( id_data ).template cast<double>() ;

      ceres::CostFunction *cost_fun =
        new ceres::AutoDiffCostFunction<RigidMotion3d3dEstimationCost, 3 , 4 , 3>( new RigidMotion3d3dEstimationCost( ref_pt, data_pt ) );
      problem.AddResidualBlock( cost_fun, nullptr , q.coeffs().data() , t.data() );
    }
  }

  ceres::Solver::Options solverOptions;
  solverOptions.minimizer_progress_to_stdout = false;
  solverOptions.logging_type                 = ceres::SILENT;
  solverOptions.linear_solver_type = ceres::DENSE_QR;

  // !!! This is important !
  // parametrization (in order to make orthogonal transformations and to use eigen ordering)
  ceres::LocalParameterization* quaternion_local_parameterization =
    new ceres::EigenQuaternionParameterization;
  problem.SetParameterization( q.coeffs().data(), quaternion_local_parameterization );

#ifdef OPENMVG_USE_OPENMP
  solverOptions.num_threads               = omp_get_max_threads();
  solverOptions.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP

  ceres::Solver::Summary summary;
  ceres::Solve( solverOptions, &problem, &summary );

  return std::make_pair( q , t );
}

/**
* @brief Computes ridig motion (R,t) : R * data + t that maps data on ref
* @param ref Reference point cloud
* @param ref_normal Reference point cloud normals
* @param data Data point cloud
* @param corresps (id of the correspondances from data to ref, -1 if no corresp)
* @return < quat , tra > Rigid motion that minimise the MSE
* use point to normal estimation
*/
template <typename Scalar>
std::pair<Eigen::Quaternion<Scalar>, Eigen::Matrix<Scalar, 3, 1>> RigidMotion3d3dEstimation<Scalar>::operator()( const point_list_type &ref,
    const point_list_type &ref_normal,
    const point_list_type &data, const std::vector<int> &corresps ) const
{
  // 1 Add residuals functions
  ceres::Problem problem;
  // Initial value (no transformation)
  Eigen::Quaternion<Scalar> q = Eigen::Quaternion<Scalar>::Identity();
  Eigen::Matrix<Scalar, 3, 1> t( 0, 0, 0 );

  for ( size_t id_data = 0; id_data < data.rows() ; ++id_data )
  {
    if ( corresps[ id_data ] >= 0 )
    {
      const int id_ref  = corresps[ id_data ];

      const Vec3 ref_pt = ref.row( id_ref ).template cast<double>() ;
      const Vec3 ref_n = ref_normal.row( id_ref ).template cast<double>() ;
      const Vec3 data_pt = data.row( id_data ).template cast<double>() ;

      ceres::CostFunction *cost_fun =
        new ceres::AutoDiffCostFunction<RigidMotion3d3dNormalEstimationCost, 1 , 4 , 3>( new RigidMotion3d3dNormalEstimationCost( ref_pt , ref_n , data_pt ) );
      problem.AddResidualBlock( cost_fun, nullptr , q.coeffs().data() , t.data() );
    }
  }

  ceres::Solver::Options solverOptions;
  solverOptions.minimizer_progress_to_stdout = false;
  solverOptions.logging_type                 = ceres::SILENT;
  solverOptions.linear_solver_type = ceres::DENSE_QR;

  // !!! This is important !
  // parametrization (in order to make orthogonal transformations and to use eigen ordering)
  ceres::LocalParameterization* quaternion_local_parameterization =
    new ceres::EigenQuaternionParameterization;
  problem.SetParameterization( q.coeffs().data(), quaternion_local_parameterization );

#ifdef OPENMVG_USE_OPENMP
  solverOptions.num_threads               = omp_get_max_threads();
  solverOptions.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP

  ceres::Solver::Summary summary;
  ceres::Solve( solverOptions, &problem, &summary );

  return std::make_pair( q , t );
}

} // namespace registration
} // namespace geometry
} // namespace openMVG

#endif