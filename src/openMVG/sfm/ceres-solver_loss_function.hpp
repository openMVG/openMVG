#ifndef CERES_SOLVER_LOSS_FUNCTION_H_
#define CERES_SOLVER_LOSS_FUNCTION_H_

#include "ceres/loss_function.h"

namespace ceres {

// from paper "A General and Adaptive Robust Loss Function", when alpha=-2
// f(x, −2, c) = 2*(x/c)^2 / [(x/c)^2+4] reproduces Geman-McClure loss
// s = x^2
class CERES_EXPORT GMCLoss : public LossFunction {
public:
  explicit GMCLoss(double a) : a_(a), b_(a * a) { }
  virtual void Evaluate(double, double*) const;

private:
  const double a_;
  // b = a^2.
  const double b_;
};

} // namespace ceres

#endif // CERES_SOLVER_LOSS_FUNCTION_H_
