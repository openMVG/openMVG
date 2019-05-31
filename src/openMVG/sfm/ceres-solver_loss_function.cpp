#include "openMVG/sfm/ceres-solver_loss_function.hpp"

namespace ceres {

void GMCLoss::Evaluate(double s, double rho[3]) const {
   const double h = s + 4.0 * b_;
   rho[0] = 2.0 * s / h;
   rho[1] = 8.0 * b_ / h / h;
   rho[2] = -16.0 * b_ / h / h / h;
}

} // namespace ceres