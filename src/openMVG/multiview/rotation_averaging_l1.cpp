// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 cDc and Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/rotation_averaging_l1.hpp"
#include "openMVG/numeric/l1_solver_admm.hpp"

#ifdef HAVE_BOOST
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
using namespace boost;
#else
#include "lemon/adaptors.h"
#include "lemon/dfs.h"
#include "lemon/kruskal.h"
#include "lemon/list_graph.h"
#include "lemon/path.h"
using namespace lemon;
#endif

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>

#include <queue>

namespace openMVG   {
namespace rotation_averaging  {
namespace l1  {

/////////////////////////

// given an array of values, compute the X84 threshold as in:
// Hampel FR, Rousseeuw PJ, Ronchetti EM, Stahel WA
// "Robust Statistics: the Approach Based on Influence Functions"
// Wiley Series in Probability and Mathematical Statistics, John Wiley & Sons, 1986
// returns the pair(median,trust_region)
// upper-bound threshold = median+trust_region
// lower-bound threshold = median-trust_region
template<typename TYPE>
inline std::pair<TYPE, TYPE>
ComputeX84Threshold(const TYPE* const values, uint32_t size, TYPE mul=TYPE(5.2))
{
  assert(size > 0);
  typename std::vector<TYPE> data(values, values+size);
  typename std::vector<TYPE>::iterator mid = data.begin() + size / 2;
  std::nth_element(data.begin(), mid, data.end());
  const TYPE median = *mid;
  // threshold = 5.2 * MEDIAN(ABS(values-median));
  for (size_t i=0; i<size; ++i)
    data[i] = std::abs(values[i]-median);
  std::nth_element(data.begin(), mid, data.end());
  return {median, mul*(*mid)};
} // ComputeX84Threshold


/////////////////////////

using Matrix3x3 = openMVG::Mat3;
using IndexArr = std::vector<uint32_t>;

// find the shortest cycle for the given graph and starting vertex
struct Node {
  using InternalType = IndexArr;
  InternalType edges; // array of vertex indices
};
using NodeArr = std::vector<Node>;

struct Link {
  uint32_t ID; // node index
  uint32_t parentID;// parent link
  inline Link(uint32_t ID_=0, uint32_t parentID_=0) : ID(ID_), parentID(parentID_) {}
};
using LinkQue = std::queue<Link>;

#ifdef HAVE_BOOST
using edge_property_t = boost::property<boost::edge_weight_t, float>;
using graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, uint32_t, edge_property_t>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using edge_iter = boost::graph_traits<graph_t>::edge_iterator;
#else
using graph_t = lemon::ListGraph;
using map_EdgeMap = graph_t::EdgeMap<double>;
#endif
using MapEdgeIJ2R = std::map<std::pair<uint32_t,uint32_t>, Matrix3x3>;

// Look for the maximum spanning tree along the graph of relative rotations
// since we look for the maximum spanning tree using a minimum spanning tree algorithm
// weight are negated.
uint32_t FindMaximumSpanningTree(const RelativeRotations& RelRs, graph_t& g, MapEdgeIJ2R& mapIJ2R, NodeArr& minGraph)
{
  assert(!RelRs.empty());
#ifdef HAVE_BOOST
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    boost::add_edge(relR.i, relR.j, - relR.weight, g);
    mapIJ2R[{relR.i, relR.j}] = relR.Rij;
    mapIJ2R[{relR.j, relR.i}] = relR.Rij.transpose();
  }
  // find the minimum spanning tree
  const size_t nViews = boost::num_vertices(g);
  minGraph.resize(nViews);
  std::vector<edge_t> spanningTree;
  boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanningTree));
  for (std::vector<edge_t>::const_iterator ei=spanningTree.begin(); ei!=spanningTree.end(); ++ei) {
    const edge_t& edge = *ei;
    minGraph[edge.m_source].edges.push_back(edge.m_target);
    minGraph[edge.m_target].edges.push_back(edge.m_source);
  }
  const size_t nEdges = spanningTree.size();
  return nEdges;
#else

  //A-- Compute the number of node we need
  std::set<uint32_t> setNodes;
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    setNodes.insert(relR.i);
    setNodes.insert(relR.j);
  }

  //B-- Create a node graph for each element of the set
  using map_NodeT = std::map<uint32_t, graph_t::Node>;
  map_NodeT map_index_to_node;
  for (const auto & iter : setNodes)
  {
    map_index_to_node[iter] = g.addNode();
  }

  //C-- Create a graph from RelRs with weighted edges
  map_EdgeMap map_edgeMap(g);
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    mapIJ2R[{relR.i, relR.j}] = relR.Rij;
    mapIJ2R[{relR.j, relR.i}] = relR.Rij.transpose();

    // add edge to the graph
    graph_t::Edge edge =  g.addEdge(map_index_to_node[relR.i], map_index_to_node[relR.j]);
    map_edgeMap[ edge ] = - relR.weight;
  }

  //D-- Compute the MST of the graph
  std::vector<graph_t::Edge> tree_edge_vec;
  lemon::kruskal(g, map_edgeMap, std::back_inserter(tree_edge_vec));

  const size_t nViews = lemon::countNodes(g);
  minGraph.resize(nViews);

  //E-- Export compute MST
  for (size_t i= 0; i < tree_edge_vec.size(); i++)
  {
    minGraph[g.id(g.u(tree_edge_vec[i]))].edges.push_back(g.id(g.v(tree_edge_vec[i])));
    minGraph[g.id(g.v(tree_edge_vec[i]))].edges.push_back(g.id(g.u(tree_edge_vec[i])));
  }
  return tree_edge_vec.size();
#endif
}
//----------------------------------------------------------------


// Filter the given relative rotations using the known global rotations
// returns the number of inliers
unsigned int FilterRelativeRotations(
  const RelativeRotations& RelRs,
  const Matrix3x3Arr& Rs,
  float threshold,
  std::vector<bool> * vec_inliers)
{
  assert(!RelRs.empty() && !Rs.empty());
  assert(threshold >= 0);
  // compute errors for each relative rotation
  std::vector<float> errors(RelRs.size());
  for (size_t r= 0; r<RelRs.size(); ++r) {
    const RelativeRotation& relR = RelRs[r];
    const Matrix3x3& Ri = Rs[relR.i];
    const Matrix3x3& Rj = Rs[relR.j];
    const Matrix3x3& Rij = relR.Rij;
    const Mat3 eRij(Rj.transpose()*Rij*Ri);
    const openMVG::Vec3 erij;
    ceres::RotationMatrixToAngleAxis((const double*)eRij.data(), (double*)erij.data());
    errors[r] = (float)erij.norm();
  }
  if (threshold == 0) {
    // estimate threshold
    const std::pair<float,float> res = ComputeX84Threshold(&errors[0], errors.size());
    threshold = res.first+res.second;
  }
  if (vec_inliers)  {
    vec_inliers->resize(RelRs.size());
  }
  // mark outliers
  unsigned int nInliers = 0;
  for (size_t r=0; r<errors.size(); ++r) {
    const bool bInlier = errors[r] < threshold;
    if (vec_inliers)
      (*vec_inliers)[r] = bInlier;
    if (bInlier)
      ++nInliers;
  }
  return nInliers;
} // FilterRelativeRotations
//----------------------------------------------------------------


double RelRotationAvgError
(
  const RelativeRotations& RelRs,
  const Matrix3x3Arr& Rs,
  double* pMin=nullptr,
  double* pMax=nullptr
)
{
#ifdef HAVE_BOOST
  boost::accumulators::accumulator_set<double,
    boost::accumulators::stats<
      boost::accumulators::tag::min,
      boost::accumulators::tag::mean,
      boost::accumulators::tag::max> > acc;

  for (int i=0; i < RelRs.size(); ++i) {
    const RelativeRotation& relR = RelRs[i];
    acc(openMVG::FrobeniusNorm(relR.Rij  - (Rs[relR.j]*Rs[relR.i].transpose())));
  }
  if (pMin)
    *pMin = boost::accumulators::min(acc);
  if (pMax)
    *pMax = boost::accumulators::max(acc);
  return boost::accumulators::mean(acc);
#else
  std::vector<double> vec_err(RelRs.size(), 0.0);
  for (size_t i=0; i < RelRs.size(); ++i) {
    const RelativeRotation& relR = RelRs[i];
    vec_err[i] = openMVG::FrobeniusNorm(relR.Rij  - (Rs[relR.j]*Rs[relR.i].transpose()));
  }
  float min, max, mean, median;
  minMaxMeanMedian(vec_err.begin(), vec_err.end(), min, max, mean, median);
  if (pMin)
    *pMin = min;
  if (pMax)
    *pMax = max;
  return mean;
#endif
}
//----------------------------------------------------------------

void InitRotationsMST
(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID
)
{
  assert(!Rs.empty());

  // -- Compute coarse global rotation estimates:
  //   - by finding the maximum spanning tree and linking the relative rotations
  //   - Initial solution is driven by relative rotations data confidence.
  graph_t g;
  MapEdgeIJ2R mapIJ2R;
  NodeArr minGraph;
  // find the Maximum Spanning Tree
  FindMaximumSpanningTree(RelRs, g, mapIJ2R, minGraph);
  g.clear();

  // start from the main view and link all views using the relative rotation estimates
  LinkQue stack;
  stack.push(Link(nMainViewID, uint32_t(0)));
  Rs[nMainViewID] = Matrix3x3::Identity();
  do {
    const Link& link = stack.front();
    const Node& node = minGraph[link.ID];

    for (Node::InternalType::const_iterator pEdge = node.edges.begin();
      pEdge != node.edges.end(); ++pEdge) {
        const size_t edge = *pEdge;
        if (edge == link.parentID) {
          // compute the global rotation for the current node
          assert(mapIJ2R.find(std::make_pair(link.parentID, link.ID)) != mapIJ2R.end());
          const Matrix3x3& Rij = mapIJ2R[{link.parentID, link.ID}];
          Rs[link.ID] = Rij * Rs[link.parentID];
        } else {
          // add edge to the processing queue
          stack.push(Link(edge, link.ID));
        }
    }
    stack.pop();
  } while(!stack.empty());
}

// Robustly estimate global rotations from relative rotations as in:
// "Efficient and Robust Large-Scale Rotation Averaging", Chatterjee and Govindu, 2013
// and detect outliers relative rotations and return them with 0 in arrInliers
bool GlobalRotationsRobust(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID,
  float threshold,
  std::vector<bool> * vec_Inliers)
{
  assert(!Rs.empty());

  // -- Compute coarse global rotation estimates:
  InitRotationsMST(RelRs, Rs, nMainViewID);

  // refine global rotations based on the relative rotations
  const bool bOk = RefineRotationsAvgL1IRLS(RelRs, Rs, nMainViewID);

  // find outlier relative rotations
  if (threshold>=0 && vec_Inliers)  {
    FilterRelativeRotations(RelRs, Rs, threshold, vec_Inliers);
  }

  return bOk;
} // GlobalRotationsRobust
//----------------------------------------------------------------

namespace internal
{

// build A in Ax=b
inline void FillMappingMatrix(
  const RelativeRotations& RelRs,
  const uint32_t nMainViewID,
  sMat& A)
{
  A.reserve(A.rows()*2); // estimate of the number of non-zeros (optional)
  sMat::Index i = 0, j = 0;
  for (size_t r=0; r<RelRs.size(); ++r) {
    const RelativeRotation& relR = RelRs[r];
    if (relR.i != nMainViewID) {
      j = 3*(relR.i<nMainViewID ? relR.i : relR.i-1);
      A.insert(i+0,j+0) = -1.0;
      A.insert(i+1,j+1) = -1.0;
      A.insert(i+2,j+2) = -1.0;
    }
    if (relR.j != nMainViewID) {
      j = 3*(relR.j<nMainViewID ? relR.j : relR.j-1);
      A.insert(i+0,j+0) = 1.0;
      A.insert(i+1,j+1) = 1.0;
      A.insert(i+2,j+2) = 1.0;
    }
    i+=3;
  }
  A.makeCompressed();
}

// compute errors for each relative rotation
inline void FillErrorMatrix(
  const RelativeRotations& RelRs,
  const Matrix3x3Arr& Rs,
  Vec & b)
{
  for (size_t r = 0; r < RelRs.size(); ++r) {
    const RelativeRotation& relR = RelRs[r];
    const Matrix3x3& Ri = Rs[relR.i];
    const Matrix3x3& Rj = Rs[relR.j];
    const Matrix3x3& Rij = relR.Rij;
    const Mat3 eRij(Rj.transpose()*Rij*Ri);
    const openMVG::Vec3 erij;
    ceres::RotationMatrixToAngleAxis((const double*)eRij.data(), (double*)erij.data());
    b.block<3,1>(3*r,0) = erij;
  }
}

// apply correction to global rotations
inline void CorrectMatrix(
  const Mat& x,
  const uint32_t nMainViewID,
  Matrix3x3Arr& Rs)
{
  for (size_t r = 0; r < Rs.size(); ++r) {
    if (r == nMainViewID)
      continue;
    Matrix3x3& Ri = Rs[r];
    const uint32_t i = (r<nMainViewID ? r : r-1);
    const openMVG::Vec3 eRid = openMVG::Vec3(x.block<3,1>(3*i,0));
    const Mat3 eRi;
    ceres::AngleAxisToRotationMatrix((const double*)eRid.data(), (double*)eRi.data());
    Ri = Ri*eRi;
  }
}

// L1RA -> L1 Rotation Averaging implementation
bool SolveL1RA
(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const sMat & A,
  const unsigned int nMainViewID
)
{
  const unsigned nObss = (unsigned)RelRs.size();
  const unsigned nVars = (unsigned)Rs.size()-1; // one view is kept constant
  const unsigned m = nObss*3;
  const unsigned n = nVars*3;

  // init x with 0 that corresponds to trusting completely the initial Ri guess
  Vec x(Vec::Zero(n)), b(m);

  // Current error and the previous one
  double e = std::numeric_limits<double>::max(), ep;
  unsigned iter = 0;
  // L1RA iterate optimization till the desired precision is reached
  do {
    // compute errors for each relative rotation
    FillErrorMatrix(RelRs, Rs, b);

    // solve the linear system using l1 norm
    L1Solver<sMat >::Options options;
    L1Solver<sMat > l1_solver(options, A);
    l1_solver.Solve(b, &x);

    ep = e; e = x.norm();
    if (ep < e)
      break;
    // apply correction to global rotations
    CorrectMatrix(x, nMainViewID, Rs);
  } while (++iter < 32 && e > 1e-5 && (ep-e)/e > 1e-2);

  std::cout << "L1RA Converged in " << iter << " iterations." << std::endl;

  return true;
}

// Iteratively Reweighted Least Squares (IRLS) implementation
bool SolveIRLS
(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const sMat & A,
  const unsigned int nMainViewID,
  const double sigma
)
{
  const unsigned nObss = (unsigned)RelRs.size();
  const unsigned nVars = (unsigned)Rs.size()-1; // one view is kept constant
  const unsigned m = nObss*3;
  const unsigned n = nVars*3;

  // init x with 0 that corresponds to trusting completely the initial Ri guess
  Vec x(Vec::Zero(n)), b(m);

  // Since the sparsity pattern will not change with each linear solve
  //  compute it once to speed up the solution time.
  using Linear_Solver_T = Eigen::SimplicialLDLT<sMat >;

  Linear_Solver_T linear_solver;
  linear_solver.analyzePattern(A.transpose() * A);
  if (linear_solver.info() != Eigen::Success) {
    std::cerr << "Cholesky decomposition failed." << std::endl;
    return false;
  }

  const double sigmaSq(Square(sigma));

  Eigen::ArrayXd errors, weights;
  Vec xp(n);
  // current error and the previous one
  double e = std::numeric_limits<double>::max(), ep;
  unsigned int iter = 0;
  do
  {
    xp = x;
    // compute errors for each relative rotation
    FillErrorMatrix(RelRs, Rs, b);

    // Compute the weights for each error term
    errors = (A * x - b).array();

    // compute robust errors using the Huber-like loss function
    weights = sigmaSq / (errors.square() + sigmaSq).square();

    // Update the factorization for the weighted values
    const sMat at_weight = A.transpose() * weights.matrix().asDiagonal();
    linear_solver.factorize(at_weight * A);
    if (linear_solver.info() != Eigen::Success) {
      std::cerr << "Failed to factorize the least squares system." << std::endl;
      return false;
    }

    // Solve the least squares problem
    x = linear_solver.solve(at_weight * b);
    if (linear_solver.info() != Eigen::Success) {
      std::cerr << "Failed to solve the least squares system." << std::endl;
      return false;
    }

    // apply correction to global rotations
    CorrectMatrix(x, nMainViewID, Rs);

    ep = e; e = (xp-x).norm();

  } while (++iter < 32 && e > 1e-5 && (ep-e)/e > 1e-2);

  std::cout << "IRLS Converged in " << iter << " iterations." << std::endl;

  return true;
}

} // namespace internal

// Refine the global rotations using to the given relative rotations, similar to:
// "Efficient and Robust Large-Scale Rotation Averaging", Chatterjee and Govindu, 2013
// L1 Rotation Averaging (L1RA) and Iteratively Reweighted Least Squares (IRLS) implementations combined
bool RefineRotationsAvgL1IRLS(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID,
  const double sigma)
{
  assert(!RelRs.empty() && !Rs.empty());
  assert(Rs[nMainViewID] == Matrix3x3::Identity());

  double fMinBefore, fMaxBefore;
  const double fMeanBefore = RelRotationAvgError(RelRs, Rs, &fMinBefore, &fMaxBefore);

  const unsigned nObss = (unsigned)RelRs.size();
  const unsigned nVars = (unsigned)Rs.size()-1; // main view is kept constant
  const unsigned m = nObss*3;
  const unsigned n = nVars*3;

  // build mapping matrix A in Ax=b
  sMat A(m, n);
  internal::FillMappingMatrix(RelRs, nMainViewID, A);

  if (!internal::SolveL1RA(RelRs, Rs, A, nMainViewID))
  {
    std::cerr << "Could not solve the L1 regression step." << std::endl;
    return false;
  }

  if (!internal::SolveIRLS(RelRs, Rs, A, nMainViewID, sigma))
  {
    std::cerr << "Could not solve the ILRS step." << std::endl;
    return false;
  }

  double fMinAfter, fMaxAfter;
  const double fMeanAfter = RelRotationAvgError(RelRs, Rs, &fMinAfter, &fMaxAfter);

  std::cout << "Refine global rotations using L1RA-IRLS and " << nObss << " relative rotations:\n"
    << " error reduced from " << fMeanBefore << "(" <<fMinBefore << " min, " << fMaxBefore << " max)\n"
    << " to " << fMeanAfter << "(" << fMinAfter << "min,"<< fMaxAfter<< "max)" << std::endl;

  return true;
} // RefineRotationsAvgL1IRLS

} // namespace l1
} // namespace rotation_averaging
} // namespace openMVG
