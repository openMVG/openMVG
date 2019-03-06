// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/clustering/kmeans.hpp"

#include "testing/testing.h"

#include <random>
#include <vector>

using namespace openMVG;
using namespace clustering;

// Define test constant:
static const int NB_CLUSTER = 3;
static const int NB_POINT = 1e4;
static const std::array<int, 3> POINTS_PER_CLUSTER =
  {NB_POINT, NB_POINT, NB_POINT};
static const std::array<KMeansInitType, 2> KMEAN_INIT_TYPES =
  {KMeansInitType::KMEANS_INIT_RANDOM, KMeansInitType::KMEANS_INIT_PP};

// Initialize NB_CLUSTER centers and POINTS_PER_CLUSTER[i] points around each centroid
// Note: Clusters and points are column based
// The function return all the points randomly initialized around the centroid.
Mat InitRandom3ClusterDataset
(
  Mat & centers,
  const int dimension
)
{
  centers.resize(dimension, NB_CLUSTER);
  centers.col(0).setConstant(-5.0);
  centers.col(1).setConstant(5.0);
  centers.col(2).setConstant(0.0);

  // Initialize a list of points around the centroids
  Mat points(dimension,
             std::accumulate(POINTS_PER_CLUSTER.cbegin(),
                             POINTS_PER_CLUSTER.cend(),
                             0));

  std::mt19937_64 rng(std::mt19937_64::default_seed);
  std::uniform_real_distribution<double> distrib(-2.0, 2.0);

  // Generate points on each centers
  int point_id = 0;
  for (int cluster_id = 0; cluster_id < NB_CLUSTER; ++cluster_id)
  {
    for (int id = 0; id < POINTS_PER_CLUSTER[cluster_id]; ++id)
    {
      Vec random_point(dimension);
      for (int i = 0; i < dimension; ++i)
        random_point(i) = distrib(rng);
      points.col(point_id) = random_point + centers.col(cluster_id);
      ++point_id;
    }
  }
  return points;
}

// Template function to convert the data_point dataset to various type
template <typename T> T ConvertMat2T(const Mat & mat);

// Predicate to test if two element are the same
bool EqualPredicate (uint32_t i, uint32_t j) {
  return (i==j);
}

// Template function to make the testing of KMeans easier
// It avoids some code duplication
template<typename ContainerType>
void KMeanTesting
(
  std::vector< uint32_t > & ids,
  std::vector< typename ContainerType::value_type > & centers,
  const int dimension,
  const KMeansInitType k_mean_init_type,
  const uint32_t k_mean_centers = 3
)
{
  // Data initialization (centers and data_points)
  Mat mat_centers;
  const Mat mat_points =
    InitRandom3ClusterDataset(mat_centers, dimension);
  const ContainerType pts =
    ConvertMat2T<ContainerType>(mat_points);

  // K-Means clustering:
  KMeans(pts, ids, centers, k_mean_centers,
         std::numeric_limits<uint32_t>::max(), k_mean_init_type );
}

// Check the result of the KMean Ids classification
#define KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids) \
{\
  /* Does the right number of cluster are found?*/\
  EXPECT_EQ(NB_CLUSTER, centers.size());\
\
  /* Does points are labelled to valid cluster ids?*/\
  {\
    std::vector< uint32_t > ids_cpy(ids.size());\
    const auto it = std::unique_copy (ids.cbegin(),ids.cend(),ids_cpy.begin(), EqualPredicate);\
    ids_cpy.resize( std::distance(ids_cpy.begin(), it) );\
    EXPECT_EQ(NB_CLUSTER, ids_cpy.size());\
  }\
\
  /* Does points are labelled to the right centroid id?*/\
  int counter = 0;\
  for (const auto nb_point_in_cluster : POINTS_PER_CLUSTER)\
  {\
    const int id_to_check = ids[counter];\
    for (int i = 0; i < nb_point_in_cluster; ++i)\
    {\
      EXPECT_EQ(id_to_check, ids[counter]);\
      ++counter;\
    }\
  }\
}

// std::vector<Vec2> specialization
template <>
std::vector< Vec2 > ConvertMat2T<std::vector< Vec2 >>
(
  const Mat & mat
)
{
  std::vector< Vec2 > pts(mat.cols());
  for (int i = 0; i < mat.cols(); ++i)
    pts[i] = mat.col(i);
  return pts;
}

TEST( clustering, threeClustersVec2 )
{
  const int dimension = 2;
  using DataPointType = Vec2;
  using ContainerType = std::vector<DataPointType>;

  for (const auto kmean_init_type : KMEAN_INIT_TYPES)
  {
    std::vector<uint32_t> ids;
    std::vector<DataPointType> centers;
    KMeanTesting<ContainerType>
    (
      ids,
      centers,
      dimension,
      kmean_init_type
    );

    KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids);

    std::cout << "Centers with kmean init type(" << (int)kmean_init_type << "): " << std::endl;
    for( const auto & it : centers )
    {
      std::cout << it.transpose() << std::endl;
    }
  }
}

// std::vector<Vec3> specialization
template <>
std::vector< Vec3 > ConvertMat2T<std::vector< Vec3 >>(const Mat & mat)
{
  std::vector< Vec3 > pts(mat.cols());
  for (int i = 0; i < mat.cols(); ++i)
    pts[i] = mat.col(i);
  return pts;
}

TEST( clustering, threeClustersVec3 )
{
  const int dimension = 3;
  using DataPointType = Vec3;
  using ContainerType = std::vector<DataPointType>;

  for (const auto kmean_init_type : KMEAN_INIT_TYPES)
  {
    std::vector<uint32_t> ids;
    std::vector<DataPointType> centers;
    KMeanTesting<ContainerType>
    (
      ids,
      centers,
      dimension,
      kmean_init_type
    );

    KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids);

    std::cout << "Centers with kmean init type(" << (int)kmean_init_type << "): " << std::endl;
    for( const auto & it : centers )
    {
      std::cout << it.transpose() << std::endl;
    }
  }
}

// std::vector<std::array<double, 5>> specialization
template <>
std::vector<std::array<double, 5>>
ConvertMat2T<std::vector<std::array<double, 5>>>(const Mat & mat)
{
  std::vector<std::array<double, 5>> pts(mat.cols());
  for (int i = 0; i < mat.cols(); ++i)
    std::copy(mat.col(i).data(), mat.col(i).data()+5, pts[i].data());
  return pts;
}

TEST( clustering, threeClustersStdArray )
{
  const int dimension = 5;
  using DataPointType = std::array<double, 5>;
  using ContainerType = std::vector<DataPointType>;

  for (const auto kmean_init_type : KMEAN_INIT_TYPES)
  {
    std::vector<uint32_t> ids;
    std::vector<DataPointType> centers;
    KMeanTesting<ContainerType>
    (
      ids,
      centers,
      dimension,
      kmean_init_type
    );

    KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids);

    std::cout << "Centers with kmean init type(" << (int)kmean_init_type << "): " << std::endl;
    for( const auto & center : centers )
    {
      for (const auto val : center)
        std::cout << val << " ";
      std::cout << std::endl;
    }
  }
}

// std::vector<std::vector<double>> specialization
template <>
std::vector<std::vector<double>>
ConvertMat2T<std::vector<std::vector<double>>>(const Mat & mat)
{
  std::vector<std::vector<double>> pts(mat.cols());
  for (int i = 0; i < mat.cols(); ++i)
    pts[i].assign(mat.col(i).data(), mat.col(i).data() + mat.rows());
  return pts;
}

TEST( clustering, threeClustersStdVector )
{
  const int dimension = 12;
  using DataPointType = std::vector<double>;
  using ContainerType = std::vector<DataPointType>;

  for (const auto kmean_init_type : KMEAN_INIT_TYPES)
  {
    std::vector<uint32_t> ids;
    std::vector<DataPointType> centers;
    KMeanTesting<ContainerType>
    (
      ids,
      centers,
      dimension,
      kmean_init_type
    );

    KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids);

    std::cout << "Centers with kmean init type(" << (int)kmean_init_type << "): " << std::endl;
    for( const auto & center : centers )
    {
      for (const auto val : center)
        std::cout << val << " ";
      std::cout << std::endl;
    }
  }
}

// std::vector<Vec> specialization
template <>
std::vector< Vec > ConvertMat2T<std::vector< Vec >>(const Mat & mat)
{
  std::vector< Vec > pts(mat.cols());
  for (int i = 0; i < mat.cols(); ++i)
    pts[i] = mat.col(i);
  return pts;
}

TEST( clustering, threeClustersEigenVec )
{
  const int dimension = 6;
  using DataPointType = Vec;
  using ContainerType = std::vector<DataPointType>;

  for (const auto kmean_init_type : KMEAN_INIT_TYPES)
  {
    std::vector<uint32_t> ids;
    std::vector<DataPointType> centers;
    KMeanTesting<ContainerType>
    (
      ids,
      centers,
      dimension,
      kmean_init_type
    );

    KMEANS_CHECK_VALIDITY(NB_CLUSTER, ids);

    std::cout << "Centers with kmean init type(" << (int)kmean_init_type << "): " << std::endl;
    for( const auto & it : centers )
    {
      std::cout << it.transpose() << std::endl;
    }
  }
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
