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

  //EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(trifocal_model_t);
  // datum_i[4 /*xy tgtx tgty*/][pp:npoints /* 3 for Chicago */]
  static void Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor);
  
  // Gabriel's comment: If bearing is the bearing vector of the camera, Vec3 should be used instead of Mat32 or use &bearing.data()[0] 
  static double Error(
    const trifocal_model_t &tt,
    const Vec &bearing_0, // x,y,tangentialx,tangentialy
    const Vec &bearing_1,
    const Vec &bearing_2,
    // --
    const Vec &pxbearing_0,  // TODO: remove. Pixel coordinates is for development/debug, we want to compute the error in pixels
    const Vec &pxbearing_1,
    const Vec &pxbearing_2,
    const double K[2][3]);
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
    : x1_(x1), x2_(x2), x3_(x3), pxx1_(pxx1), pxx2_(pxx2), pxx3_(pxx3), K_(K) {}

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
    // std::cerr << "\n" << sample << "\n";
    return ErrorArg::Error(model, 
        x1_.col(sample), x2_.col(sample), x3_.col(sample), 
        pxx1_.col(sample), pxx2_.col(sample), pxx3_.col(sample), K_);
  }

  /// Number of putative point
  size_t NumSamples() const { return static_cast<size_t>(x1_.cols()); }

  /// Compute a model on sampled datum_
  static void Solve(const Mat &x1, const Mat &x2, const Mat &x3, vector<Model> *models) {
    Solver::Solve(x1, x2, x3, models); // By offering this, Kernel types can be passed to templates.
  }
  
protected:
  const Mat &x1_, &x2_, &x3_; // corresponding point of the trifical configuration
  // x_i[4 /*xy tgtx tgty*/][npts /* total number of tracks */]
  // -- Below in pixels for debug
  const Mat &pxx1_, &pxx2_, &pxx3_;
  const double (*K_)[3]; // pointer to 2x3 array
};


} // namespace trifocal3pt
#endif // trifocal_sample_h_