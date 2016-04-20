
// Copyright (c) 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//------------------
//-- Bibliography --
//------------------
//- [1] "Global Multiple-View Color Consistency."
//- Authors: Pierre MOULON, Bruno DUISIT and Pascal MONASSE.
//- Date: November 2013.
//- Conference: CVMP.

//-- Linear programming
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#include "openMVG/linearProgramming/bisectionLP.hpp"

namespace openMVG {
namespace lInfinity {

struct relativeColorHistogramEdge
{
  size_t I,J;
  std::vector<size_t> histoI, histoJ;

  relativeColorHistogramEdge() = default ;

  relativeColorHistogramEdge(
    size_t i, size_t j,
    const std::vector<size_t> & histogramI,
    const std::vector<size_t> & histogramJ):
      I(i), J(j),
      histoI(histogramI),
      histoJ(histogramJ)
    { }
};

namespace histogram {

// Normalize a distribution function
template<typename T>
static void normalizeHisto(const std::vector<T> & vec_df, std::vector<double> & vec_normalized_df)
{
  double totalCount = static_cast<double>(std::accumulate(vec_df.begin(), vec_df.end(), 0));
  vec_normalized_df.resize(vec_df.size(), 0.0);
  for(size_t i=0; i<vec_df.size(); i++)
    vec_normalized_df[i] = vec_df[i] / totalCount;
}

// Compute cumulative distribution functions (cdf)
template<typename T>
static void cdf(const std::vector<T> & vec_df, std::vector<T> & vec_cdf)
{
  vec_cdf = vec_df;
  for(size_t i=1; i<vec_cdf.size(); i++)
      vec_cdf[i] = vec_cdf[i] + vec_cdf[i-1];
}

}; // namespace histogram

// Implementation of the formula (1) of [1] with 10 quantiles.
//-- L_infinity alignment of pair of histograms over a graph thanks to a linear program.
static void Encode_histo_relation(
    const size_t nImage,
    const std::vector<relativeColorHistogramEdge > & vec_relativeHistograms,
    const std::vector<size_t> & vec_indexToFix,
    sRMat & A, Vec & C,
    std::vector<linearProgramming::LP_Constraints::eLP_SIGN> & vec_sign,
    std::vector<double> & vec_costs,
    std::vector< std::pair<double,double> > & vec_bounds)
{
  const size_t Nima = (size_t) nImage;
  const size_t Nrelative = vec_relativeHistograms.size();

# define GVAR(i) (2*(i))
# define OFFSETVAR(i) (2*(i)+1)
# define GAMMAVAR (2*Nima)

  const size_t nbQuantile = 10;

  const size_t Nconstraint = nbQuantile * Nrelative * 2;
  const size_t NVar = 2 * Nima+ 1;

  A.resize(Nconstraint, NVar);

  C.resize(Nconstraint, 1);
  C.fill(0.0);
  vec_sign.resize(Nconstraint);

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(NVar);
  fill( vec_bounds.begin(), vec_bounds.end(),
    std::make_pair((double)-1e+30, (double)1e+30));

  // Set gain as positive values
  for (size_t i = 0; i < Nima; ++i)
  {
    vec_bounds[GVAR(i)].first = 0.0;
  }

  //-- Fix the required image to known gain and offset value
  for (std::vector<size_t>::const_iterator iter = vec_indexToFix.begin();
    iter != vec_indexToFix.end(); ++iter)
  {
    vec_bounds[GVAR(*iter)] = std::make_pair(1.0, 1.0);      // gain = 1.0
    vec_bounds[OFFSETVAR(*iter)] = std::make_pair(0.0, 0.0); // offset = 0
  }

  // Setup gamma >= 0
  vec_bounds[GAMMAVAR].first = 0.0;

  //-- Minimize gamma
  vec_costs.resize(NVar);
  std::fill(vec_costs.begin(), vec_costs.end(), 0.0);
  vec_costs[GAMMAVAR] = 1.0;
  //--

  size_t rowPos = 0;
  double incrementPourcentile = 1./(double) nbQuantile;

  for (size_t i = 0; i < Nrelative; ++i)
  {
    std::vector<relativeColorHistogramEdge>::const_iterator iter = vec_relativeHistograms.begin();
    std::advance(iter, i);

    const relativeColorHistogramEdge & edge = *iter;

    //-- compute the two cumulated and normalized histogram

    const std::vector< size_t > & vec_histoI = edge.histoI;
    const std::vector< size_t > & vec_histoJ = edge.histoJ;

    const size_t nBuckets = vec_histoI.size();

    // Normalize histogram
    std::vector<double> ndf_I(nBuckets), ndf_J(nBuckets);
    histogram::normalizeHisto(vec_histoI, ndf_I);
    histogram::normalizeHisto(vec_histoJ, ndf_J);

    // Compute cumulative distribution functions (cdf)
    std::vector<double> cdf_I(nBuckets), cdf_J(nBuckets);
    histogram::cdf(ndf_I, cdf_I);
    histogram::cdf(ndf_J, cdf_J);

    double currentPourcentile = 5./100.;

    //-- Compute pourcentile and their positions
    std::vector<double> vec_pourcentilePositionI, vec_pourcentilePositionJ;
    vec_pourcentilePositionI.reserve(1.0/incrementPourcentile);
    vec_pourcentilePositionJ.reserve(1.0/incrementPourcentile);

    std::vector<double>::const_iterator cdf_I_IterBegin = cdf_I.begin();
    std::vector<double>::const_iterator cdf_J_IterBegin = cdf_J.begin();
    while( currentPourcentile < 1.0)
    {
      std::vector<double>::const_iterator iterFI = std::lower_bound(cdf_I.begin(), cdf_I.end(), currentPourcentile);
      const size_t positionI = std::distance(cdf_I_IterBegin, iterFI);

      std::vector<double>::const_iterator iterFJ = std::lower_bound(cdf_J.begin(), cdf_J.end(), currentPourcentile);
      const size_t positionJ = std::distance(cdf_J_IterBegin, iterFJ);

      vec_pourcentilePositionI.push_back(positionI);
      vec_pourcentilePositionJ.push_back(positionJ);

      currentPourcentile += incrementPourcentile;
    }

    //-- Add the constraints:
    // pos * ga + offa - pos * gb - offb <= gamma
    // pos * ga + offa - pos * gb - offb >= - gamma

    for(size_t k = 0; k < vec_pourcentilePositionI.size(); ++k)
    {
      A.coeffRef(rowPos, GVAR(edge.I)) = vec_pourcentilePositionI[k];
      A.coeffRef(rowPos, OFFSETVAR(edge.I)) = 1.0;

      A.coeffRef(rowPos, GVAR(edge.J)) = - vec_pourcentilePositionJ[k];
      A.coeffRef(rowPos, OFFSETVAR(edge.J)) = - 1.0;

      // - gamma (side change)
      A.coeffRef(rowPos, GAMMAVAR) = -1;
      // <= gamma
      vec_sign[rowPos] = linearProgramming::LP_Constraints::LP_LESS_OR_EQUAL;
      C(rowPos) = 0;
      ++rowPos;

      A.coeffRef(rowPos, GVAR(edge.I)) = vec_pourcentilePositionI[k];
      A.coeffRef(rowPos, OFFSETVAR(edge.I)) = 1.0;

      A.coeffRef(rowPos, GVAR(edge.J)) = - vec_pourcentilePositionJ[k];
      A.coeffRef(rowPos, OFFSETVAR(edge.J)) = - 1.0;

      // + gamma (side change)
      A.coeffRef(rowPos, GAMMAVAR) = 1;
      // >= - gamma
      vec_sign[rowPos] = linearProgramming::LP_Constraints::LP_GREATER_OR_EQUAL;
      C(rowPos) = 0;
      ++rowPos;
    }
  }
#undef GVAR
#undef OFFSETVAR
#undef GAMMAVAR
}

struct ConstraintBuilder_GainOffset
{
  ConstraintBuilder_GainOffset(
    const std::vector<relativeColorHistogramEdge > & vec_relativeHistograms,
    const std::vector<size_t> & vec_indexToFix):
    _vec_relative(vec_relativeHistograms),
    _vec_indexToFix(vec_indexToFix)
  {
    //Count the number of images
    std::set<size_t> countSet;
    for (int i = 0; i  < _vec_relative.size(); ++i)
    {
      countSet.insert(_vec_relative[i].I);
      countSet.insert(_vec_relative[i].J);
    }
    _Nima = countSet.size();
  }

  /// Setup constraints for the translation and structure problem,
  ///  in the LP_Constraints object.
  bool Build(linearProgramming::LP_Constraints_Sparse & constraint)
  {
    Encode_histo_relation(
      _Nima,
      _vec_relative,
      _vec_indexToFix,
      constraint.constraint_mat_,
      constraint.constraint_objective_,
      constraint.vec_sign_,
      constraint.vec_cost_,
      constraint.vec_bounds_);

    // it's a minimization problem over the gamma variable
    constraint.bminimize_ = true;

    //-- Setup additional information about the Linear Program constraint
    constraint.nbParams_ = _Nima * 2 + 1;
    return true;
  }
  // Internal data
  size_t _Nima;
  const std::vector< relativeColorHistogramEdge > & _vec_relative;
  const std::vector<size_t> & _vec_indexToFix;
};


}; // namespace lInfinity
}; // namespace openMVG
