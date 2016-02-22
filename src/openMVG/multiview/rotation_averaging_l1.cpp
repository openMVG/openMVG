
// Copyright (c) 2014 cDc and Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/rotation_averaging_l1.hpp"

#ifdef HAVE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/foreach.hpp>
using namespace boost;
#else
#include "lemon/list_graph.h"
#include "lemon/kruskal.h"
#include "lemon/adaptors.h"
#include "lemon/dfs.h"
#include "lemon/path.h"
using namespace lemon;
#endif


#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include <map>
#include <queue>
#include <stdint.h>

namespace openMVG   {
namespace rotation_averaging  {
namespace l1  {

// Minimum l1 error approximation:
//
// Let A be a M x N matrix with full rank. Given y of R^M, the problem
// (PA) minx ||y - Ax||1
// finds the vector x of R^N such that the error y - Ax has minimum l1 norm
// (i.e. we are asking that the difference between Ax and y be sparse).
// This problem arises in the context of channel coding
// (see E. J. Candes and T. Tao. "Decoding by linear programming". in IEEE Trans.
// Inform. Theory, December 2005).
//
// Suppose we have a channel code that produces a codeword c = Ax for a message x. The
// message travels over the channel, and has an unknown number of its entries corrupted.
// The decoder observes y = c + e, where e is the corruption. If e is sparse enough, then
// the decoder can use (PA) to recover x exactly. When x, A, y have real-valued entries,
// (PA) can be recast as an LP.
//
template<typename MATRIX_TYPE>
inline bool TRobustRegressionL1PD(
  const MATRIX_TYPE& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& y,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& xp,
  REAL pdtol, unsigned pdmaxiter)
{
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> Vector;
  const unsigned M = (unsigned)y.size();
  const unsigned N = (unsigned)xp.size();
  assert(A.rows() == M && A.cols() == N);

  const REAL alpha(0.01);
  const REAL beta(0.5);
  const REAL mu(10);

  Vector x(xp);
  Vector Ax(A*x);
  Vector tmpM1(y-Ax);
  Vector tmpM2(-tmpM1);
  Vector tmpM3(tmpM1.cwiseAbs()), tmpM4(M);
  Vector u = (tmpM3*REAL(0.95)).array() + tmpM3.maxCoeff()*REAL(0.10);
  Vector fu1 = tmpM2-u;
  Vector fu2 = tmpM1-u;

  Vector lamu1(M), lamu2(M);
  for (unsigned i=0; i<M; ++i) {
    lamu1(i) = -1.0/fu1(i);
    lamu2(i) = -1.0/fu2(i);
  }
  const MATRIX_TYPE At(A.transpose());
  Vector Atv(At*(lamu1-lamu2));
  REAL AtvNormSq = Atv.squaredNorm();
  Vector rdual((-lamu1-lamu2).array() + REAL(1));
  REAL rdualNormSq = rdual.squaredNorm();

  Vector w2(M), sig1(M), sig2(M), sigx(M), dx(N), up(N), Atdv(N);
  Vector Axp(M), Atvp(M);
  Vector &Adx(sigx), &du(w2), &w1p(dx);
  Matrix H11p(N,N);
  Vector &dlamu1(tmpM3), &dlamu2(tmpM4);
  for (unsigned pditer=0; pditer<pdmaxiter; ++pditer) {
    // surrogate duality gap
    const REAL sdg(-(fu1.dot(lamu1) + fu2.dot(lamu2)));
    if (sdg < pdtol)
      break;
    const REAL tau(mu*2*M/sdg);
    const REAL inv_tau = REAL(-1)/tau;
    tmpM1 = (-lamu1.cwiseProduct(fu1)).array() + inv_tau;
    tmpM2 = (-lamu2.cwiseProduct(fu2)).array() + inv_tau;
    const REAL resnorm = sqrt(AtvNormSq + rdualNormSq + tmpM1.squaredNorm() + tmpM2.squaredNorm());

    for (unsigned i=0; i<M; ++i) {
      REAL& tmpM3i = tmpM3(i);
      tmpM3i = inv_tau/fu1(i);
      REAL& tmpM4i = tmpM4(i);
      tmpM4i = inv_tau/fu2(i);
      w2(i) = tmpM3i + tmpM4i - REAL(1);
    }

    tmpM1 = lamu1.cwiseQuotient(fu1);
    tmpM2 = lamu2.cwiseQuotient(fu2);
    sig1 = -tmpM1 - tmpM2;
    sig2 = tmpM1 - tmpM2;
    sigx = sig1 - sig2.cwiseAbs2().cwiseQuotient(sig1);

    H11p = At*(Eigen::DiagonalMatrix<REAL,Eigen::Dynamic>(sigx)*A);
    w1p = At*(tmpM4 - tmpM3 - (sig2.cwiseQuotient(sig1).cwiseProduct(w2)));

    // optimized solver as A is positive definite and symmetric
    dx = H11p.ldlt().solve(w1p);

    Adx = A*dx;

    du = (w2 - sig2.cwiseProduct(Adx)).cwiseQuotient(sig1);

    dlamu1 = -tmpM1.cwiseProduct(Adx-du) - lamu1 + tmpM3;
    dlamu2 =  tmpM2.cwiseProduct(Adx+du) - lamu2 + tmpM4;
    Atdv = At*(dlamu1-dlamu2);

    // make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
    REAL s(1);
    for (unsigned i=0; i<M; ++i) {
      REAL& dlamu1i = dlamu1(i);
      if (dlamu1i < 0) {
        const REAL tmp = -lamu1(i)/dlamu1i;
        if (s > tmp)
          s = tmp;
      }
      REAL& dlamu2i = dlamu2(i);
      if (dlamu2i < 0) {
        const REAL tmp = -lamu2(i)/dlamu2i;
        if (s > tmp)
          s = tmp;
      }
    }
    for (unsigned i=0; i<M; ++i) {
      REAL& Adxi = Adx(i);
      REAL& dui = du(i);
      const REAL Adx_du = Adxi-dui;
      if (Adx_du > 0) {
        const REAL tmp = -fu1(i)/Adx_du;
        if (s > tmp)
          s = tmp;
      }
      const REAL _Adx_du = -Adxi-dui;
      if (_Adx_du > 0) {
        const REAL tmp = -fu2(i)/_Adx_du;
        if (s > tmp)
          s = tmp;
      }
    }
    s *= REAL(0.99);

    // backtrack
    lamu1 += s*dlamu1;  lamu2 += s*dlamu2;
    rdual = (-lamu1-lamu2).array() + REAL(1);
    rdualNormSq = rdual.squaredNorm();
    bool suffdec = false;
    unsigned backiter = 0;
    do {
      xp = x + s*dx;  up = u + s*du;
      Axp = Ax + s*Adx;  Atvp = Atv + s*Atdv;
      fu1 = Axp - y - up;  fu2 = -Axp + y - up;
      AtvNormSq = Atvp.squaredNorm();
      tmpM1 = (-lamu1.cwiseProduct(fu1)).array() + inv_tau;
      tmpM2 = (-lamu2.cwiseProduct(fu2)).array() + inv_tau;
      const REAL newresnorm = sqrt(AtvNormSq + rdualNormSq + tmpM1.squaredNorm() + tmpM2.squaredNorm());
      suffdec = (newresnorm <= (REAL(1)-alpha*s)*resnorm);
      s = beta*s;
      if (++backiter > 32) {
        //("error: stuck backtracking, returning last iterate"); // see Section 4 of notes for more information
        xp.swap(x);
        return false;
      }
    } while (!suffdec);

    // next iteration
    x.swap(xp);  u.swap(up);
    Ax.swap(Axp);  Atv.swap(Atvp);
  }
  return true;
}

bool RobustRegressionL1PD(
  const Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& b,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& x,
  REAL pdtol, unsigned pdmaxiter)
{
  return TRobustRegressionL1PD(A, b, x, pdtol, pdmaxiter);
}
bool RobustRegressionL1PD(
  const Eigen::SparseMatrix<REAL, Eigen::ColMajor>& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& b,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& x,
  REAL pdtol, unsigned pdmaxiter)
{
  return TRobustRegressionL1PD(A, b, x, pdtol, pdmaxiter);
}

/*----------------------------------------------------------------*/

// Iteratively Re-weighted Least Squares (IRLS) implementation
template<typename MATRIX_TYPE>
inline bool TIterativelyReweightedLeastSquares(
  const MATRIX_TYPE& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& b,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& x,
  REAL sigma, REAL eps)
{
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> Vector;
  const unsigned m = (unsigned)b.size();
  const unsigned n = (unsigned)x.size();
  assert(A.rows() == m && A.cols() == n);

  // iterate optimization till the desired precision is reached
  Vector xp(n), e(m);
  const REAL sigmaSq(Square(sigma));
  unsigned iter = 0;
  REAL delta = std::numeric_limits<REAL>::max(), deltap;
  do {
    xp = x;
    // compute error vector
    e = A*x-b;
    // compute robust errors using the Huber-like loss function
    for (unsigned i=0; i<m; ++i) {
      REAL& err = e(i);
      const REAL errSq(Square(err));
      err = sigmaSq / (errSq + sigmaSq);
    }
    // solve the linear system using l2 norm
    const MATRIX_TYPE AtF(A.transpose()*e.asDiagonal());
    const Eigen::LDLT<Matrix> solver(AtF*A); // compute the Cholesky decomposition
    if (solver.info() != Eigen::Success) {
      std::cerr << "error: decomposing linear system failed" << std::endl;
      return false;
    }
    x = solver.solve(AtF*b);
    if (solver.info() != Eigen::Success) {
      std::cerr << "error: solving linear system failed" << std::endl;
      return false;
    }
    if (++iter > 32)
      break;
    deltap = delta; delta = (xp-x).norm();
  } while (delta > eps && (deltap-delta)/delta > 1e-2);
  return true;
}

bool IterativelyReweightedLeastSquares(
  const Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& b,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& x,
  REAL sigma, REAL eps)
{
  return TIterativelyReweightedLeastSquares(A, b, x, sigma, eps);
}
bool IterativelyReweightedLeastSquares(
  const Eigen::SparseMatrix<REAL, Eigen::ColMajor>& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& b,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& x,
  REAL sigma, REAL eps)
{
  return TIterativelyReweightedLeastSquares(A, b, x, sigma, eps);
}

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
ComputeX84Threshold(const TYPE* const values, size_t size, TYPE mul=TYPE(5.2))
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
  return std::make_pair(median, mul*(*mid));
} // ComputeX84Threshold


/////////////////////////

typedef openMVG::Mat3 Matrix3x3;
typedef std::vector<size_t> IndexArr;

// find the shortest cycle for the given graph and starting vertex
struct Node {
  typedef IndexArr InternalType;
  InternalType edges; // array of vertex indices
};
typedef std::vector<Node> NodeArr;

struct Link {
  size_t ID; // node index
  size_t parentID;// parent link
  inline Link(size_t _ID=0, size_t _parentID=0) : ID(_ID), parentID(_parentID) {}
};
typedef std::queue<Link> LinkQue;

#ifdef HAVE_BOOST
typedef boost::property<boost::edge_weight_t, float> edge_property_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, size_t, edge_property_t> graph_t;
typedef graph_t::vertex_descriptor vertex_t;
typedef graph_t::edge_descriptor edge_t;
typedef boost::graph_traits<graph_t>::edge_iterator edge_iter;
#else
typedef lemon::ListGraph graph_t;
typedef graph_t::EdgeMap<double> map_EdgeMap;
#endif
typedef std::map<std::pair<size_t,size_t>, Matrix3x3> MapEdgeIJ2R;

// Look for the maximum spanning tree along the graph of relative rotations
// since we look for the maximum spanning tree using a minimum spanning tree algorithm
// weight are negated.
size_t FindMaximumSpanningTree(const RelativeRotations& RelRs, graph_t& g, MapEdgeIJ2R& mapIJ2R, NodeArr& minGraph)
{
  assert(!RelRs.empty());
#ifdef HAVE_BOOST
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    boost::add_edge(relR.i, relR.j, - relR.weight, g);
    mapIJ2R[std::make_pair(relR.i, relR.j)] = relR.Rij;
    mapIJ2R[std::make_pair(relR.j, relR.i)] = relR.Rij.transpose();
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
  std::set<size_t> setNodes;
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    setNodes.insert(relR.i);
    setNodes.insert(relR.j);
  }

  //B-- Create a node graph for each element of the set
  typedef std::map<size_t, graph_t::Node> map_Size_t_Node;
  map_Size_t_Node map_size_t_to_node;
  for (std::set<size_t>::const_iterator iter = setNodes.begin();
    iter != setNodes.end();
    ++iter)
  {
    map_size_t_to_node[*iter] = g.addNode();
  }

  //C-- Create a graph from RelRs with weighted edges
  map_EdgeMap map_edgeMap(g);
  for (size_t p = 0; p < RelRs.size(); ++p) {
    const RelativeRotation& relR = RelRs[p];
    mapIJ2R[std::make_pair(relR.i, relR.j)] = relR.Rij;
    mapIJ2R[std::make_pair(relR.j, relR.i)] = relR.Rij.transpose();

    // add edge to the graph
    graph_t::Edge edge =  g.addEdge(map_size_t_to_node[relR.i], map_size_t_to_node[relR.j]);
    map_edgeMap[ edge ] = - relR.weight;
  }

  //D-- Compute the MST of the graph
  std::vector<graph_t::Edge> tree_edge_vec;
  lemon::kruskal(g, map_edgeMap, std::back_inserter(tree_edge_vec));

  const size_t nViews = lemon::countNodes(g);
  minGraph.resize(nViews);

  //E-- Export compute MST
  for(size_t i= 0 ; i < tree_edge_vec.size(); i++)
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
  for(int r= 0; r<RelRs.size(); ++r) {
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
  for(int r=0; r<errors.size(); ++r) {
    bool bInlier = (errors[r] < threshold);
    if (vec_inliers)
      (*vec_inliers)[r] = bInlier;
    if (bInlier)
      ++nInliers;
  }
  return nInliers;
} // FilterRelativeRotations
//----------------------------------------------------------------


REAL RelRotationAvgError(const RelativeRotations& RelRs, const Matrix3x3Arr& Rs, REAL* pMin=NULL, REAL* pMax=NULL)
{
#ifdef HAVE_BOOST
  boost::accumulators::accumulator_set<REAL,
    boost::accumulators::stats<
      boost::accumulators::tag::min,
      boost::accumulators::tag::mean,
      boost::accumulators::tag::max> > acc;

  for(int i=0; i < RelRs.size(); ++i) {
    const RelativeRotation& relR = RelRs[i];
    acc(openMVG::FrobeniusNorm(relR.Rij  - (Rs[relR.j]*Rs[relR.i].transpose())));
  }
  if (pMin)
    *pMin = boost::accumulators::min(acc);
  if (pMax)
    *pMax = boost::accumulators::max(acc);
  return boost::accumulators::mean(acc);
#else
  std::vector<REAL> vec_err(RelRs.size(), REAL(0.0));
  for(int i=0; i < RelRs.size(); ++i) {
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
  const size_t nMainViewID
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
  stack.push(Link(nMainViewID, size_t(0)));
  Rs[nMainViewID] = Matrix3x3::Identity();
  do {
    const Link& link = stack.front();
    const Node& node = minGraph[link.ID];

    for(Node::InternalType::const_iterator pEdge = node.edges.begin();
      pEdge != node.edges.end(); ++pEdge) {
        const size_t edge = *pEdge;
        if (edge == link.parentID) {
          // compute the global rotation for the current node
          assert(mapIJ2R.find(std::make_pair(link.parentID, link.ID)) != mapIJ2R.end());
          const Matrix3x3& Rij = mapIJ2R[std::make_pair(link.parentID, link.ID)];
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
  const size_t nMainViewID,
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


// build A in Ax=b
inline void _FillMappingMatrix(
  const RelativeRotations& RelRs,
  const size_t nMainViewID,
  Eigen::SparseMatrix<REAL,Eigen::ColMajor>& A)
{
  A.reserve(A.rows()*2); // estimate of the number of non-zeros (optional)
  Eigen::SparseMatrix<REAL,Eigen::ColMajor>::Index i = 0, j = 0;
  for(int r=0; r<RelRs.size(); ++r) {
    const RelativeRotation& relR = RelRs[r];
    if (relR.i != nMainViewID) {
      j = 3*(relR.i<nMainViewID ? relR.i : relR.i-1);
      A.insert(i+0,j+0) = REAL(-1);
      A.insert(i+1,j+1) = REAL(-1);
      A.insert(i+2,j+2) = REAL(-1);
    }
    if (relR.j != nMainViewID) {
      j = 3*(relR.j<nMainViewID ? relR.j : relR.j-1);
      A.insert(i+0,j+0) = REAL(1);
      A.insert(i+1,j+1) = REAL(1);
      A.insert(i+2,j+2) = REAL(1);
    }
    i+=3;
  }
  A.makeCompressed();
}

// compute errors for each relative rotation
inline void _FillErrorMatrix(
  const RelativeRotations& RelRs,
  const Matrix3x3Arr& Rs,
  Eigen::Matrix<REAL,Eigen::Dynamic,1>& b)
{
  for (size_t r = 0; r < RelRs.size(); ++r) {
    const RelativeRotation& relR = RelRs[r];
    const Matrix3x3& Ri = Rs[relR.i];
    const Matrix3x3& Rj = Rs[relR.j];
    const Matrix3x3& Rij = relR.Rij;
    const Mat3 eRij(Rj.transpose()*Rij*Ri);
    const openMVG::Vec3 erij;
    ceres::RotationMatrixToAngleAxis((const double*)eRij.data(), (double*)erij.data());
    b.block<3,1>(3*r,0) = openMVG::Vec3(erij*relR.weight);
  }
}

// apply correction to global rotations
inline void _CorrectMatrix(
  const Eigen::Matrix<REAL,Eigen::Dynamic,1>& x,
  const size_t nMainViewID,
  Matrix3x3Arr& Rs)
{
  for (size_t r = 0; r < Rs.size(); ++r) {
    if (r == nMainViewID)
      continue;
    Matrix3x3& Ri = Rs[r];
    const size_t i = (r<nMainViewID ? r : r-1);
    openMVG::Vec3 eRid = openMVG::Vec3(x.block<3,1>(3*i,0));
    const Mat3 eRi;
    ceres::AngleAxisToRotationMatrix((const double*)eRid.data(), (double*)eRi.data());
    Ri = Ri*eRi;
  }
}

// Refine the global rotations using to the given relative rotations, similar to:
// "Efficient and Robust Large-Scale Rotation Averaging", Chatterjee and Govindu, 2013
// L1 Rotation Averaging (L1RA) and Iteratively Reweighted Least Squares (IRLS) implementations combined
bool RefineRotationsAvgL1IRLS(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const size_t nMainViewID,
  REAL sigma)
{
  assert(!RelRs.empty() && !Rs.empty());
  assert(Rs[nMainViewID] == Matrix3x3::Identity());

  REAL fMinBefore, fMaxBefore, fMeanBefore = RelRotationAvgError(RelRs, Rs, &fMinBefore, &fMaxBefore);

  const unsigned nObss = (unsigned)RelRs.size();
  const unsigned nVars = (unsigned)Rs.size()-1; // main view is kept constant
  const unsigned m = nObss*3;
  const unsigned n = nVars*3;

  // build mapping matrix A in Ax=b
  Eigen::SparseMatrix<REAL,Eigen::ColMajor> A(m, n);
  _FillMappingMatrix(RelRs, nMainViewID, A);

  // init x with 0 that corresponds to trusting completely the initial Ri guess
  Vec x(Vec::Zero(n)), b(m);

  // L1RA iterate optimization till the desired precision is reached
  REAL e = std::numeric_limits<REAL>::max(), ep;
  unsigned iter1 = 0;
  do {
    // compute errors for each relative rotation
    _FillErrorMatrix(RelRs, Rs, b);
    // solve the linear system using l1 norm
    if (!RobustRegressionL1PD(A, b, x)) {
      std::cerr << "error: l1 robust regression failed." << std::endl;
      return false;
    }
    ep = e; e = x.norm();
    if (ep < e)
      break;
    // apply correction to global rotations
    _CorrectMatrix(x, nMainViewID, Rs);
  } while (++iter1 < 32 && e > 1e-5 && (ep-e)/e > 1e-2);
  // IRLS iterate optimization till the desired precision is reached
  x.setZero();
  e = std::numeric_limits<REAL>::max();
  unsigned iter2 = 0;
  do {
    // compute errors for each relative rotation
    _FillErrorMatrix(RelRs, Rs, b);
    // solve the linear system using l2 norm
    if (!IterativelyReweightedLeastSquares(A, b, x, sigma)) {
      std::cerr << "error: l2 iterative regression failed" << std::endl;
      return false;
    }
    ep = e; e = x.norm();
    if (ep < e)
      break;
    // apply correction to global rotations
    _CorrectMatrix(x, nMainViewID, Rs);
  } while (++iter2 < 32 && e > 1e-5 && (ep-e)/e > 1e-2);

  REAL fMinAfter, fMaxAfter, fMeanAfter = RelRotationAvgError(RelRs, Rs, &fMinAfter, &fMaxAfter);

  std::cout << "Refine global rotations using L1RA-IRLS and " << nObss << " relative rotations:\n"
    << " error reduced from " << fMeanBefore << "(" <<fMinBefore << " min, " << fMaxBefore << " max)\n"
    << " to " << fMeanAfter << "(" << fMinAfter << "min,"<< fMaxAfter<< "max)\n"
    << " in " << iter1 << "+" << iter2 << "=" << iter1+iter2 << " iterations" << std::endl;

  return true;
} // RefineRotationsAvgL1IRLS

} // namespace l1
} // namespace rotation_averaging
} // namespace openMVG

