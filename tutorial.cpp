#include "ceres/ceres.h"
#include "glog/logging.h"

#include <vector>

// class QuadraticCostFunction : public ceres::SizedCostFunction<1, 1> {
//  public:
//   virtual ~QuadraticCostFunction() {}
//   virtual bool Evaluate(double const* const* parameters,
//                         double* residuals,
//                         double** jacobians) const {
//     const double x = parameters[0][0];
//     residuals[0] = 10 - x;

//     // Compute the Jacobian if asked for.
//     if (jacobians != nullptr && jacobians[0] != nullptr) {
//       jacobians[0][0] = -1;
//     }
//     return true;
//   }
// };

struct F1 {
  template <typename T>
  bool operator()(const T* const x1, const T* const x2, T* residual) const {
    // f1 = x1 + 10 * x2;
    residual[0] = x1[0] + 10.0 * x2[0];
    return true;
  }
};
struct F2 {
  template <typename T>
  bool operator()(const T* const x3, const T* const x4, T* residual) const {
    // f2 = sqrt(5) (x3 - x4)
    residual[0] = sqrt(5.0) * (x3[0] - x4[0]);
    return true;
  }
};
struct F3 {
  template <typename T>
  bool operator()(const T* const x2, const T* const x3, T* residual) const {
    // f3 = (x2 - 2 x3)^2
    residual[0] = (x2[0] - 2.0 * x3[0]) * (x2[0] - 2.0 * x3[0]);
    return true;
  }
};
struct F4 {
  template <typename T>
  bool operator()(const T* const x1, const T* const x4, T* residual) const {
    // f4 = sqrt(10) (x1 - x4)^2
    residual[0] = sqrt(10.0) * (x1[0] - x4[0]) * (x1[0] - x4[0]);
    return true;
  }
};

DEFINE_string(minimizer,
              "trust_region",
              "Minimizer type to use, choices are: line_search & trust_region");

int main(int argc, char* argv[])
{
    GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // The variable to solve for with its initial value.
    double x1 = 3.0;
    double x2 = -1.0;
    double x3 = 0.0;
    double x4 = 1.0;

    // Build the problem.
    ceres::Problem problem;
    // Add residual terms to the problem using the autodiff
    // wrapper to get the derivatives automatically. The parameters, x1 through
    // x4, are modified in place.
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<F1, 1, 1, 1>(new F1), nullptr, &x1, &x2);
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<F2, 1, 1, 1>(new F2), nullptr, &x3, &x4);
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<F3, 1, 1, 1>(new F3), nullptr, &x2, &x3);
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<F4, 1, 1, 1>(new F4), nullptr, &x1, &x4);

 

  ceres::Solver::Options options;
  LOG_IF(FATAL,
         !ceres::StringToMinimizerType(CERES_GET_FLAG(FLAGS_minimizer),
                                       &options.minimizer_type))
      << "Invalid minimizer: " << CERES_GET_FLAG(FLAGS_minimizer)
      << ", valid options are: trust_region and line_search.";
  options.max_num_iterations = 100;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;
  // clang-format off
  std::cout << "Initial x1 = " << x1
            << ", x2 = " << x2
            << ", x3 = " << x3
            << ", x4 = " << x4
            << "\n";
  // clang-format on
  // Run the solver!
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << "\n";
  // clang-format off
  std::cout << "Final x1 = " << x1
            << ", x2 = " << x2
            << ", x3 = " << x3
            << ", x4 = " << x4
            << "\n";
  // clang-format on








    // // Set up the only cost function (also known as residual). This uses
    // // auto-differentiation to obtain the derivative (jacobian).
    // // ceres::CostFunction* cost_function =
    // //     new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);

    // ceres::CostFunction* cost_function = new QuadraticCostFunction;
    // problem.AddResidualBlock(cost_function, nullptr, &x);

    // // Run the solver!
    // ceres::Solver::Options options;
    // options.linear_solver_type = ceres::DENSE_QR;
    // options.minimizer_progress_to_stdout = true;
    // ceres::Solver::Summary summary;
    // ceres::Solve(options, &problem, &summary);

    // std::cout << summary.FullReport() << "\n";
    // std::cout << "x : " << initial_x
    //           << " -> " << x << "\n";
  return 0;
}