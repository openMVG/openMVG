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

    // Compute residual projection cost given two points
    struct RigidMotion3d3dEstimationCost
    {
    public:
      RigidMotion3d3dEstimationCost( const Vec3 &ref, const Vec3 &data )
          : m_ref( ref ),
            m_data( data )
      {
      }

      template <typename T>
      bool operator()( const T *const x, T *residual ) const
      {
        // x : (R[0;8] | t[0;2])
        Eigen::Matrix<T, 3, 1> src;
        src << T( m_data[ 0 ] ), T( m_data[ 1 ] ), T( m_data[ 2 ] );

        const T tx = x[ 0 ] * src[ 0 ] + x[ 1 ] * src[ 1 ] + x[ 2 ] * src[ 2 ] + x[ 9 ];
        const T ty = x[ 3 ] * src[ 0 ] + x[ 4 ] * src[ 1 ] + x[ 5 ] * src[ 2 ] + x[ 10 ];
        const T tz = x[ 6 ] * src[ 0 ] + x[ 7 ] * src[ 1 ] + x[ 8 ] * src[ 2 ] + x[ 11 ];

        Eigen::Matrix<T, 3, 1> tra;
        tra << tx, ty, tz;

        // 2 - compute distance between transformed and
        Eigen::Matrix<T, 3, 1> ref;
        ref << T( m_ref[ 0 ] ), T( m_ref[ 1 ] ), T( m_ref[ 2 ] );

        residual[ 0 ] = ref[ 0 ] - tra[ 0 ];
        residual[ 1 ] = ref[ 1 ] - tra[ 1 ];
        residual[ 2 ] = ref[ 2 ] - tra[ 2 ];

        return true;
      }

    private:
      Vec3 m_ref;
      Vec3 m_data;
    };

    template <typename Scalar>
    class RigidMotion3d3dEstimation
    {
      using point_list_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>;

    public:
      RigidMotion3d3dEstimation();

      /**
      * @brief Computes ridig motion (R,t) : R * data + t that maps data on ref 
      * @param ref Reference point cloud 
      * @param data Data point cloud 
      * @return R,t Rigid motion that minimise the MSE   
      */
      std::pair<Mat3, Vec3> operator()( const point_list_type &ref, const point_list_type &data );

      /**
      * @brief Computes ridig motion (R,t) : R * data + t that maps data on ref 
      * @param ref Reference point cloud 
      * @param data Data point cloud 
      * @param corresps (id of the correspondances from data to ref, -1 if no corresp)
      * @return R,t Rigid motion that minimise the MSE   
      */
      std::pair<Mat3, Vec3> operator()( const point_list_type &ref, const point_list_type &data, const std::vector<int> &corresps );

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
      * @return R,t Rigid motion that minimise the MSE   
      */
    template <typename Scalar>
    std::pair<Mat3, Vec3> RigidMotion3d3dEstimation<Scalar>::operator()( const RigidMotion3d3dEstimation<Scalar>::point_list_type &ref, const RigidMotion3d3dEstimation<Scalar>::point_list_type &data )
    {
      // 1 Add residuals functions
      ceres::Problem problem;
      // Initial value (no transformation)
      Scalar x[ 12 ];
      x[ 0 ] = 1.0;
      x[ 1 ] = 0.0;
      x[ 2 ] = 0.0;

      x[ 3 ] = 0.0;
      x[ 4 ] = 1.0;
      x[ 5 ] = 0.0;

      x[ 6 ] = 0.0;
      x[ 7 ] = 0.0;
      x[ 8 ] = 1.0;

      x[ 9 ]  = 0.0;
      x[ 10 ] = 0.0;
      x[ 11 ] = 0.0;

      for ( int id_obs = 0; id_obs < ref.rows(); ++id_obs )
      {
        Vec3 ref_pt( ref( id_obs, 0 ), ref( id_obs, 1 ), ref( id_obs, 2 ) );
        Vec3 data_pt( data( id_obs, 0 ), data( id_obs, 1 ), data( id_obs, 2 ) );

        ceres::CostFunction *cost_fun =
            new ceres::AutoDiffCostFunction<RigidMotion3d3dEstimationCost, 3, 12>( new RigidMotion3d3dEstimationCost( ref_pt, data_pt ) );
        problem.AddResidualBlock( cost_fun, NULL, x );
      }

      ceres::Solver::Options solverOptions;
      solverOptions.minimizer_progress_to_stdout = false;
      solverOptions.logging_type                 = ceres::SILENT;
      // Since the problem is sparse, use a sparse solver iff available
      if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::SUITE_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::CX_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::EIGEN_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else
      {
        solverOptions.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
      }
#ifdef OPENMVG_USE_OPENMP
      solverOptions.num_threads               = omp_get_max_threads();
      solverOptions.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP

      ceres::Solver::Summary summary;
      ceres::Solve( solverOptions, &problem, &summary );

      Mat3 r;

      r( 0, 0 ) = x[ 0 ];
      r( 0, 1 ) = x[ 1 ];
      r( 0, 2 ) = x[ 2 ];

      r( 1, 0 ) = x[ 3 ];
      r( 1, 1 ) = x[ 4 ];
      r( 1, 2 ) = x[ 5 ];

      r( 2, 0 ) = x[ 6 ];
      r( 2, 1 ) = x[ 7 ];
      r( 2, 2 ) = x[ 8 ];

      Vec3 t;
      t[ 0 ] = x[ 9 ];
      t[ 1 ] = x[ 10 ];
      t[ 2 ] = x[ 11 ];

      return std::make_pair( r, t );
    }

    /**
      * @brief Computes ridig motion (R,t) : R * data + t that maps data on ref 
      * @param ref Reference point cloud 
      * @param data Data point cloud 
      * @param corresps (id of the correspondances from data to ref, -1 if no corresp)
      * @return R,t Rigid motion that minimise the MSE   
      */
    template <typename Scalar>
    std::pair<Mat3, Vec3> RigidMotion3d3dEstimation<Scalar>::operator()( const point_list_type &ref,
                                                                         const point_list_type &data,
                                                                         const std::vector<int> &corresps )
    {
      // 1 Add residuals functions
      ceres::Problem problem;
      // Initial value (no transformation)
      Scalar x[ 12 ];
      x[ 0 ] = 1.0;
      x[ 1 ] = 0.0;
      x[ 2 ] = 0.0;

      x[ 3 ] = 0.0;
      x[ 4 ] = 1.0;
      x[ 5 ] = 0.0;

      x[ 6 ] = 0.0;
      x[ 7 ] = 0.0;
      x[ 8 ] = 1.0;

      x[ 9 ]  = 0.0;
      x[ 10 ] = 0.0;
      x[ 11 ] = 0.0;

      for ( size_t id_corresp = 0; id_corresp < corresps.size(); ++id_corresp )
      {
        if ( corresps[ id_corresp ] >= 0 )
        {
          const int id_data = id_corresp;
          const int id_ref  = corresps[ id_corresp ];

          Vec3 ref_pt( ref( id_ref, 0 ), ref( id_ref, 1 ), ref( id_ref, 2 ) );
          Vec3 data_pt( data( id_data, 0 ), data( id_data, 1 ), data( id_data, 2 ) );

          ceres::CostFunction *cost_fun =
              new ceres::AutoDiffCostFunction<RigidMotion3d3dEstimationCost, 3, 12>( new RigidMotion3d3dEstimationCost( ref_pt, data_pt ) );
          problem.AddResidualBlock( cost_fun, NULL, x );
        }
      }

      ceres::Solver::Options solverOptions;
      solverOptions.minimizer_progress_to_stdout = false;
      solverOptions.logging_type                 = ceres::SILENT;
      // Since the problem is sparse, use a sparse solver iff available
      if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::SUITE_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::CX_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else if ( ceres::IsSparseLinearAlgebraLibraryTypeAvailable( ceres::EIGEN_SPARSE ) )
      {
        solverOptions.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
        solverOptions.linear_solver_type                 = ceres::SPARSE_NORMAL_CHOLESKY;
      }
      else
      {
        solverOptions.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
      }
#ifdef OPENMVG_USE_OPENMP
      solverOptions.num_threads               = omp_get_max_threads();
      solverOptions.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP

      ceres::Solver::Summary summary;
      ceres::Solve( solverOptions, &problem, &summary );

      Mat3 r;

      r( 0, 0 ) = x[ 0 ];
      r( 0, 1 ) = x[ 1 ];
      r( 0, 2 ) = x[ 2 ];

      r( 1, 0 ) = x[ 3 ];
      r( 1, 1 ) = x[ 4 ];
      r( 1, 2 ) = x[ 5 ];

      r( 2, 0 ) = x[ 6 ];
      r( 2, 1 ) = x[ 7 ];
      r( 2, 2 ) = x[ 8 ];

      Vec3 t;
      t[ 0 ] = x[ 9 ];
      t[ 1 ] = x[ 10 ];
      t[ 2 ] = x[ 11 ];

      return std::make_pair( r, t );
    }

  } // namespace registration
} // namespace geometry
} // namespace openMVG

#endif