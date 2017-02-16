#ifndef OPENMVG_GEOMETRY_REGISTRATION_ICP_HPP
#define OPENMVG_GEOMETRY_REGISTRATION_ICP_HPP

#include "openMVG/geometry/kd_tree_3d.hpp"
#include "openMVG/geometry/registration/rigid_motion_3d_3d_estimation.hpp"
#include "openMVG/numeric/numeric.h"

#include <flann/flann.hpp>

#include <memory>

namespace openMVG
{
namespace geometry
{
  namespace registration
  {

    /**
    * @brief Compute centroid of a point cloud 
    * @TODO : look for an existing function or an eigen shortcut 
    */
    template <typename Scalar>
    Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> PCCentroid( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data )
    {
      Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> center;
      center.fill( Scalar( 0 ) );
      for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
      {
        center[ 0 ] += data( id_pt, 0 );
        center[ 1 ] += data( id_pt, 1 );
        center[ 2 ] += data( id_pt, 2 );
      }

      return center / static_cast<Scalar>( data.rows() );
    }

    // compute centers of target and data (restricted to pairs of bc)
    template <typename Scalar>
    std::pair<Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor>, Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor>> PCCentroid(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data,
        const std::vector<int> &biunique_corresp )
    {
      Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> tgt_center;
      Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> dtd_center;

      tgt_center.fill( Scalar( 0 ) );
      dtd_center.fill( Scalar( 0 ) );

      int nb_valid = 0;
      for ( int id_pt = 0; id_pt < biunique_corresp.size(); ++id_pt )
      {
        const int id_data   = id_pt;
        const int id_target = biunique_corresp[ id_data ];
        if ( id_target >= 0 )
        {
          tgt_center[ 0 ] += target( id_target, 0 );
          tgt_center[ 1 ] += target( id_target, 1 );
          tgt_center[ 2 ] += target( id_target, 2 );

          dtd_center[ 0 ] += data( id_data, 0 );
          dtd_center[ 1 ] += data( id_data, 1 );
          dtd_center[ 2 ] += data( id_data, 2 );
          ++nb_valid;
        }
      }

      tgt_center /= static_cast<Scalar>( nb_valid );
      dtd_center /= static_cast<Scalar>( nb_valid );

      return std::make_pair( tgt_center, dtd_center );
    }

    /**
    * @brief Compute distance between two points 
    * @TODO : look for an existing function or an eigen shortcut 
    */
    template <typename Scalar>
    Scalar PDistance( const Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> &p1, const Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> &p2 )
    {
      const Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> d = p1 - p2;

      return d[ 0 ] * d[ 0 ] +
             d[ 1 ] * d[ 1 ] +
             d[ 2 ] * d[ 2 ];
    }

    /**
    * @brief Compute threshold value as defined in BC-ICP method 
    */
    template <typename Scalar>
    Scalar Threshold( const int Nmc,
                      const Scalar lambda,
                      const Scalar lambda_c,
                      const Scalar s,
                      const Scalar c,
                      const Scalar meanSD )
    {
      if ( lambda > lambda_c )
      {
        return std::pow( Nmc, lambda ) * meanSD + s * c * c;
      }
      else
      {
        return meanSD;
      }
    }

    /**
    * @brief Transform data set using a specified transformation 
    */
    template <typename Scalar>
    void Transform( Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data, const Mat4 &tra )
    {
      for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
      {
        const Scalar x = data( id_pt, 0 );
        const Scalar y = data( id_pt, 1 );
        const Scalar z = data( id_pt, 2 );

        Vec3 pt_tra( x * tra( 0, 0 ) + y * tra( 0, 1 ) + z * tra( 0, 2 ) + tra( 0, 3 ),
                     x * tra( 1, 0 ) + y * tra( 1, 1 ) + z * tra( 1, 2 ) + tra( 1, 3 ),
                     x * tra( 2, 0 ) + y * tra( 2, 1 ) + z * tra( 2, 2 ) + tra( 2, 3 ) );
        //        const Scalar w = x * tra( 3, 0 ) + y * tra( 3, 1 ) + z * tra( 3, 2 ) + tra( 3, 3 );

        data( id_pt, 0 ) = pt_tra[ 0 ]; // / w;
        data( id_pt, 1 ) = pt_tra[ 1 ]; // / w;
        data( id_pt, 2 ) = pt_tra[ 2 ]; // / w;
      }
    }

    template <typename Scalar>
    void ComputeMinimumDistance( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
                                 const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data,
                                 std::vector<int> &index,
                                 std::vector<Scalar> &distance )
    {
      index.resize( data.rows() );
      distance.resize( data.rows() );

      for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
      {
        Scalar minimum_distance = std::numeric_limits<Scalar>::max();
        int id                  = -1;

        for ( int id_ref = 0; id_ref < target.rows(); ++id_ref )
        {
          const Scalar delta[] = {target( id_ref, 0 ) - data( id_pt, 0 ),
                                  target( id_ref, 1 ) - data( id_pt, 1 ),
                                  target( id_ref, 2 ) - data( id_pt, 2 )};
          const Scalar dist = delta[ 0 ] * delta[ 0 ] + delta[ 1 ] * delta[ 1 ] + delta[ 2 ] * delta[ 2 ];
          if( dist < minimum_distance )
          {
            minimum_distance = dist ; 
            id = id_ref ; 
          }
        }

        distance[ id_pt ] = minimum_distance;
        index[ id_pt ]    = id;
      }
    }

    /**
  * @brief Given two sets of points: target and data 
  * This function computes rigid transformation that maps model on target minimizing MSE distance 
  * @param target Target shape : a matrix of 3d points (one point per row)
  * @param data data shape : a matrix of 3d points (one point per row)
  * @param nb_iteration Maximum number of iteration 
  * @param mse_threshold Threshold use to stop computation 
  * @param[out] t Translation transformation (3d vector) 
  * @param[out] R Rotation transformation (3x3 rotation matrix)
  * @note This is an implementation based on the BC-ICP algorithm 
  */
    template <typename Scalar>
    void ICP( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
              const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data_,
              const unsigned long max_nb_iteration,
              const Scalar mse_threshold,
              openMVG::Vec3 &t,
              openMVG::Mat3 &R )
    {

      // Build Kd-Tree
      KDTree3d<Scalar> tree( target );

      unsigned long id_iteration = 0;
      Scalar cur_mse             = std::numeric_limits<Scalar>::max();
      int Nmc                    = 10 ;

      // Working sample
      Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> data = data_;

      Mat4 final_tra;
      final_tra.setIdentity();

      std::vector<bool> already_used( target.rows() );
      std::vector<int> biunique_corresp( data.rows() );
      std::vector<Scalar> biunique_distance( data.rows() );

      Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices;
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dists;

      double old_inliner_ratio = 0.0;

      while ( id_iteration < max_nb_iteration && cur_mse > mse_threshold )
      {
        std::fill( already_used.begin(), already_used.end(), false );
        std::fill( biunique_corresp.begin(), biunique_corresp.end(), -1 );
        std::fill( biunique_distance.begin(), biunique_distance.end(), Scalar(-1.0) );

        // Find Nmc nearest neighbors
        tree.Search( data, Nmc, indices, dists );

        // Compute binunique correspondance for each points
        for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
        {
          for ( int id_n = 0; id_n < Nmc; ++id_n )
          {
            const int id_point = indices( id_pt, id_n );

            if ( ( id_point >= 0 ) && ( id_point < target.rows() ) && ( !already_used[ id_point ] ) )
            {
              // We've found one !
              biunique_distance[ id_pt ] = dists( id_pt, id_n );
              biunique_corresp[ id_pt ]  = id_point;
              already_used[ id_point ]   = true;
              break;
            }
          }
        }

        // Mean squared distance (and nb inliers)
        Scalar meanSD = Scalar( 0 );
        int NbIn      = 0;
        for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
        {
          if ( biunique_corresp[ id_pt ] >= 0 )
          {
            meanSD += biunique_distance[ id_pt ];
            ++NbIn;
          }
        }
        meanSD /= static_cast<Scalar>( NbIn );

        const std::pair<Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor>,
                        Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor>>
            centers = PCCentroid( target, data, biunique_corresp );

        const Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> pc1 = centers.first;
        const Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> pc2 = centers.second;

        const Scalar c        = PDistance( pc1, pc2 );
        const Scalar s        = 1.0;
        const int nb_outliers = data.rows() - NbIn;
        const Scalar lambda   = static_cast<Scalar>( nb_outliers ) / static_cast<Scalar>( data.rows() );
        const Scalar lambda_c = 0.1;
        const Scalar t        = Threshold( Nmc, lambda, lambda_c, s, c, meanSD );

        // Filter points based on threshold
        int nb_in_prev = NbIn ; 
        for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
        {
          if ( biunique_corresp[ id_pt ] >= 0 )
          {
            if ( biunique_distance[ id_pt ] > t )
            {
              biunique_corresp[ id_pt ]  = -1;
              biunique_distance[ id_pt ] = -1.0;
            }
            else 
            {
              ++nb_in_prev ; 
            }
          }
        }

        // Compute rigid transformation using only pairs of biunique
        // 1 - Build residual fonctor
        // 2 - Run optimisation
        RigidMotion3d3dEstimation<Scalar> estimator;
        // Compute MSE error
        std::pair<Mat3, Vec3> tra = estimator( target, data, biunique_corresp );
        Mat4 tmp;
        tmp( 0, 0 ) = tra.first( 0, 0 );
        tmp( 0, 1 ) = tra.first( 0, 1 );
        tmp( 0, 2 ) = tra.first( 0, 2 );
        tmp( 0, 3 ) = tra.second[ 0 ];

        tmp( 1, 0 ) = tra.first( 1, 0 );
        tmp( 1, 1 ) = tra.first( 1, 1 );
        tmp( 1, 2 ) = tra.first( 1, 2 );
        tmp( 1, 3 ) = tra.second[ 1 ];

        tmp( 2, 0 ) = tra.first( 2, 0 );
        tmp( 2, 1 ) = tra.first( 2, 1 );
        tmp( 2, 2 ) = tra.first( 2, 2 );
        tmp( 2, 3 ) = tra.second[ 2 ];

        tmp( 3, 0 ) = 0.0;
        tmp( 3, 1 ) = 0.0;
        tmp( 3, 2 ) = 0.0;
        tmp( 3, 3 ) = 1.0;

        // Update data points
        Transform( data, tmp );
        // Compute MSE between pairs of biunique
        cur_mse      = 0.0;
        int nb_valid = 0;
        for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
        {
          const int corresp = biunique_corresp[ id_pt ];
          if ( corresp >= 0 )
          {
            const Scalar d[] = {data( id_pt, 0 ) - target( corresp, 0 ),
                                data( id_pt, 1 ) - target( corresp, 1 ),
                                data( id_pt, 2 ) - target( corresp, 2 )};

            // Add distance between the two points
            cur_mse += d[ 0 ] * d[ 0 ] + d[ 1 ] * d[ 1 ] + d[ 2 ] * d[ 2 ];
            ++nb_valid;
          }
        }
        cur_mse /= static_cast<Scalar>( nb_valid );
        std::cout << "Mse after transformation : " << cur_mse << std::endl;

        // Update final transformation and data point set
        final_tra = tmp * final_tra;

        const double cur_inlier_ratio = NbIn / data.rows();
        const double delta_ratio      = cur_inlier_ratio - old_inliner_ratio;
        if ( delta_ratio > 0.1 && Nmc > 1 )
        {
          --Nmc;
        }

        ++id_iteration;
      }

      // Compute final transformation
      R( 0, 0 ) = final_tra( 0, 0 );
      R( 0, 1 ) = final_tra( 0, 1 );
      R( 0, 2 ) = final_tra( 0, 2 );

      R( 1, 0 ) = final_tra( 1, 0 );
      R( 1, 1 ) = final_tra( 1, 1 );
      R( 1, 2 ) = final_tra( 1, 2 );

      R( 2, 0 ) = final_tra( 2, 0 );
      R( 2, 1 ) = final_tra( 2, 1 );
      R( 2, 2 ) = final_tra( 2, 2 );

      t[ 0 ] = final_tra( 0, 3 );
      t[ 1 ] = final_tra( 1, 3 );
      t[ 2 ] = final_tra( 2, 3 );
    }

  } // namespace registration
} // namespace geometry
} // namespace openMVG

#endif