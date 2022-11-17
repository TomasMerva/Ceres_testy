#include "ceres/ceres.h"
#include "glog/logging.h"
#include <iomanip>
#include <vector>
#include <array>

// class KOMO_FeatureCost : public ceres::CostFunction {
//     public:
//         KOMO_FeatureCost() {
//             set_num_residuals(1);
//             *mutable_parameter_block_sizes() = std::vector<int32_t>(1, 3);
//         }
//         bool Evaluate(double const* const* parameters,
//                               double* residuals,
//                               double** jacobians) const {
//             std::cout << parameters[0][0] << std::endl;
//             residuals[0] = 10 - parameters[0][0];

//             // Compute the Jacobian if asked for.
//             if (jacobians != nullptr && jacobians[0] != nullptr && jacobians[1] != nullptr && jacobians[2] != nullptr) {
//                 // jacobians[id_parameterblocku][id_res*block_size + id_var]
//                 // res[0]
//                 memset(jacobians[0], 0, 3*sizeof(jacobians[0]));
//                 jacobians[0][0] = -1.0;
//                 // jacobians[0][1] = 0.0;
//                 // jacobians[0][2] = 0.0;
//             }
//             return true;
//         }
// };

struct CostFunctor {
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        std::cout << x[0] << " " << x[1] << "\n";
        residual[0] = 10.0 - x[0];
        return true;
    }
};

int main(int argc, char* argv[])
{
    // std::vector<double *> x;
    // x.push_back( new double[2]{1, 2} );
    // x.push_back( new double[2]{3, 4} );

    // Eigen::ArrayXd x_old(4);
    // memcpy(x_old.data(), *x.data(), 4*sizeof(double));

    // std::cout << Eigen::VectorXd(x_old).transpose() << "\n";
    std::vector<std::vector<double>> x;
    x.push_back(std::vector<double>{1, 2});
    x.push_back(std::vector<double>{3, 4});

    ceres::Problem problem;
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 1, 2>(new CostFunctor); 
    problem.AddResidualBlock(cost_function, nullptr, x[0].data());

    double cost = 0.0;
    problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    std::cout << "evaluated cost = " << cost << std::endl;
















    // google::InitGoogleLogging(argv[0]);
    // uint num_timesteps = 20;
    // uint num_timestep_var = 3;

    // // generate q0 linspace
    // double xt[3] = {1, 4, 9};
    // double xt1[3] = {12, 45, 19};
    // double xt2[3] = {2, 12, 73};


    // // create vector of pointers to arrays representing timesteps
    // std::vector<double *> x{xt, xt1, xt2};
    // ceres::Problem problem;

    // ceres::CostFunction* cost_function = new KOMO_FeatureCost;
    // problem.AddResidualBlock(cost_function, nullptr, xt);
    // // problem.AddResidualBlock(new ceres::AutoDiffCostFunction<KOMO, 3, 3, 3, 3>(new KOMO), nullptr, xt, xt1, xt2);

    // double cost = 0.0;
    // problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

    //  // Run the solver!
    // ceres::Solver::Options options;
    // options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    // options.minimizer_progress_to_stdout = true;
    // options.max_num_iterations = 1;
    // options.parameter_tolerance = 1e-4;
    // ceres::Solver::Summary summary;
    // ceres::Solve(options, &problem, &summary);
    // std::cout << summary.FullReport() << "\n";

    return 0;
}

///////////////////////////////////////////////////////

