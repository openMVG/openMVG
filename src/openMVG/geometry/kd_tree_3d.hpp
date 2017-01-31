#ifndef OPENMVG_GEOMETRY_KDTREE3D_HPP
#define OPENMVG_GEOMETRY_KDTREE3D_HPP

#include "openMVG/numeric/numeric.h"

#include <flann/flann.hpp>

#include <memory>

namespace openMVG
{
namespace geometry
{
  /**
   * @brief Efficient data structure for searching neighbors in a 3d point cloud 
   * @TODO : find a better class name 
   * @tparam Scalar Scalar type used to define points coordinates 
   * @tparam Metric Metric used to compute distance between points 
   */
  template <typename Scalar, typename Metric = flann::L2<Scalar>>
  struct KDTree3d
  {
  public:
    using DistanceType = typename Metric::ResultType;

    /**
      * @brief Ctr 
      * @param set Set of point used for searching 
      */
    KDTree3d( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &set );

    /**
      * @brief Destructor 
      */
    ~KDTree3d();

    /**
     * @brief Search for the nearest neighbors of a given point
     * @param query Query point 
     * @param[out] index Index of the nearest point 
     * @param[out] dist Distance between query and nearest point 
     * @retval true if search is done without any problem 
     * @retval false if there was an error suring search 
     */
    bool Search( const Vec3 &query, int &index, DistanceType &dist ) const;

    /**
     * @brief Search for the nearest neighbor of a list of points 
     * @param query list of query points (one per row) 
     * @param[out] indices Index of the nearest point 
     * @param[out] dist Distance between query and it's nearest point  
     * @retval true if search is done without any problem 
     * @retval false if there was an error suring search 
     */
    bool Search( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &query,
                 std::vector<int> &indices,
                 std::vector<DistanceType> &dists ) const;

    /**
     * @brief Search for the N nearest neighbors of a given point 
     * @param query Query point
     * @param Nb number of point to search  
     * @param[out] indices Index of the nearest points 
     * @param[out] dist Distance between query and nearest points 
     * @retval true if search is done without any problem 
     * @retval false if there was an error suring search 
     */
    bool Search( const Vec3 &query, const int Nb, std::vector<int> &indices, std::vector<DistanceType> &dist ) const;

    /**
    * @brief Search for the N nearest neighbors of a list of point 
    * @param query List of query point (one per row) 
    * @param Nb Number of neighbors to search 
    * @param[out] indices of the nearest points (One row per query point, Nb columns)
    * @param[out] dists of the nearest points (One row per query point, Nb columns)
    * @retval true if search is done without any problem 
    * @retval false if there was an error suring search 
    */
    bool Search( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &query,
                 const int Nb,
                 Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &indices,
                 Eigen::Matrix<DistanceType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &dists ) const;

  private:
    /// The dataset in flann format
    std::unique_ptr<flann::Matrix<Scalar>> target_set_;
    /// The kdtree in flann format
    std::unique_ptr<flann::Index<Metric>> target_index_;
  };

  /**
   * @brief Efficient data structure for searching neighbors in a 3d point cloud 
   * @TODO : find a better class name 
   * @tparam Scalar Scalar type used to define points coordinates 
   */
  template <typename Scalar, typename Metric>
  KDTree3d<Scalar, Metric>::KDTree3d( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &set )
  {
    // Build the index
    target_set_.reset( new flann::Matrix<Scalar>( (Scalar *)set.data(), set.rows(), 3 ) );
    // Use 100% precision (an exact Kd Tree)
    target_index_.reset( new flann::Index<Metric>( *target_set_, flann::AutotunedIndexParams( 1.0 ) ) );
    target_index_->buildIndex();
  }

  /**
   * @brief Destructor 
   */
  template <typename Scalar, typename Metric>
  KDTree3d<Scalar, Metric>::~KDTree3d()
  {
    target_set_.reset();
    target_index_.reset();
  }

  /**
   * @brief Search for the nearest neighbors
   * @param query Query point 
   * @param[out] index Index of the nearest point 
   * @param[out] dist Distance between query and nearest point 
   */
  template <typename Scalar, typename Metric>
  bool KDTree3d<Scalar, Metric>::Search( const Vec3 &query, int &index, DistanceType &dist ) const
  {
    // Prepare data
    int *indexPtr   = &index;
    Scalar *distPtr = &dist;
    flann::Matrix<Scalar> queries( const_cast<Scalar *>( query.data() ), 1, 3 );
    flann::Matrix<int> indices( indexPtr, 1, 1 );
    flann::Matrix<DistanceType> dists( distPtr, 1, 1 );

    // Perform search
    return target_index_->knnSearch( queries, indices, dists, 1, flann::SearchParams( flann::FLANN_CHECKS_AUTOTUNED ) ) > 0;
  }

  /**
     * @brief Search for the nearest neighbor of a list of points 
     * @param query list of query points (one per row) 
     * @param[out] indices Index of the nearest point 
     * @param[out] dist Distance between query and it's nearest point  
     * @retval true if search is done without any problem 
     * @retval false if there was an error suring search 
     */
  template <typename Scalar, typename Metric>
  bool KDTree3d<Scalar, Metric>::Search( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &query,
                                         std::vector<int> &indices,
                                         std::vector<DistanceType> &dists ) const
  {
    // Prepare data
    indices.resize( query.rows() );
    dists.resize( query.rows() );

    int *indexPtr   = &( indices[ 0 ] );
    Scalar *distPtr = &( dists[ 0 ] );
    flann::Matrix<Scalar> queries( const_cast<Scalar *>( query.data() ), query.rows(), 3 );
    flann::Matrix<int> findices( indexPtr, query.rows(), 1 );
    flann::Matrix<DistanceType> fdists( distPtr, query.rows(), 1 );

    // Perform search
    return target_index_->knnSearch( queries, findices, fdists, 1, flann::SearchParams( flann::FLANN_CHECKS_AUTOTUNED ) ) > 0;
  }

  /**
   * @brief Search for the N nearest neighbors 
   * @param query Query point
   * @param Nb number of point to search  
   * @param[out] indices Index of the nearest points 
   * @param[out] dist Distance between query and nearest points 
   */
  template <typename Scalar, typename Metric>
  bool KDTree3d<Scalar, Metric>::Search( const Vec3 &query, const int Nb, std::vector<int> &indices, std::vector<DistanceType> &dist ) const
  {
    indices.resize( Nb );
    dist.resize( Nb );

    // Perform search
    // Prepare data
    int *indexPtr   = &( indices[ 0 ] );
    Scalar *distPtr = &( dist[ 0 ] );
    flann::Matrix<Scalar> queries( const_cast<Scalar *>( query.data() ), 1, 3 );
    flann::Matrix<int> findices( indexPtr, 1, Nb );
    flann::Matrix<DistanceType> dists( distPtr, 1, Nb );

    // Perform search
    return target_index_->knnSearch( queries, findices, dists, Nb, flann::SearchParams( flann::FLANN_CHECKS_AUTOTUNED ) ) > 0;
  }

  /**
    * @brief Search for the N nearest neighbors of a list of point 
    * @param query List of query point (one per row) 
    * @param Nb Number of neighbors to search 
    * @param[out] indices of the nearest points (One row per query point, Nb columns)
    * @param[out] dists of the nearest points (One row per query point, Nb columns)
    * @retval true if search is done without any problem 
    * @retval false if there was an error suring search 
    */
  template <typename Scalar, typename Metric>
  bool KDTree3d<Scalar, Metric>::Search( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &query,
                                         const int Nb,
                                         Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &indices,
                                         Eigen::Matrix<DistanceType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &dists ) const
  {
    indices.resize( query.rows(), Nb );
    dists.resize( query.rows(), Nb );

    // Prepare data
    int *indexPtr   = indices.data();
    Scalar *distPtr = dists.data();
    flann::Matrix<Scalar> queries( const_cast<Scalar *>( query.data() ), query.rows(), 3 );
    flann::Matrix<int> findices( indexPtr, query.rows(), Nb );
    flann::Matrix<DistanceType> fdists( distPtr, query.rows(), Nb );

    // Perform search
    return target_index_->knnSearch( queries, findices, fdists, Nb, flann::SearchParams( flann::FLANN_CHECKS_AUTOTUNED ) ) > 0;
  }

} // namespace geometry
} // namespace openMVG

#endif