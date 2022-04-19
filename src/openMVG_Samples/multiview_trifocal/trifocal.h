//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 11:55:58 -03 2021
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
#ifndef trifocal_solver_h_
#define trifocal_solver_h_

#include <iostream>
#include <array>
#include <vector>
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"


EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(std::array<openMVG::Mat34,3>)

namespace trifocal3pt {
  
using namespace std;
using namespace openMVG;

//------------------------------------------------------------------------------
struct Trifocal3PointPositionTangentialSolver {
  using trifocal_model_t = std::array<Mat34, 3>;
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };

  static void Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor);
  
  static double Error(
    const trifocal_model_t &tt,
    const Vec &bearing_0, // x,y,tangentialx,tangentialy
    const Vec &bearing_1,
    const Vec &bearing_2);
};

//------------------------------------------------------------------------------
template<typename SolverArg,
         typename ErrorArg,
         typename ModelArg = Trifocal3PointPositionTangentialSolver::trifocal_model_t>
class ThreeViewKernel {
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;
  /// The minimal number of point required for the model estimation
  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  /// The number of models that the minimal solver could return.
  enum { MAX_MODELS = Solver::MAX_MODELS };
  
  ThreeViewKernel(
      const Mat &x1, const Mat &x2, const Mat &x3, 
      const Mat &pxx1, const Mat &pxx2, const Mat &pxx3, const double K[2][3]) 
    : x1_(x1), x2_(x2), x3_(x3) {}

  /// Extract required sample and fit model(s) to the sample
  void Fit(const vector<uint32_t> &samples, vector<Model> *models) const {
    const auto
      x1 = ExtractColumns(x1_, samples),
      x2 = ExtractColumns(x2_, samples),
      x3 = ExtractColumns(x3_, samples);
    Solver::Solve(x1, x2, x3, models);
  }
  
  /// Return the error associated to the model and sample^nth point
  double Error(uint32_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1_.col(sample), x2_.col(sample), x3_.col(sample));
  }

  /// Number of putative point
  size_t NumSamples() const { return static_cast<size_t>(x1_.cols()); }

  /// Compute a model on sampled datum_
  static void Solve(const Mat &x1, const Mat &x2, const Mat &x3, vector<Model> *models) {
    Solver::Solve(x1, x2, x3, models); // By offering this, Kernel types can be passed to templates.
  }
  
protected:
  const Mat &x1_, &x2_, &x3_; // corresponding point of the trifocal configuration
};
} // namespace trifocal3pt
#endif // trifocal_sample_h_
