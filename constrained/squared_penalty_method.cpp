#include "ceres/ceres.h"
#include "glog/logging.h"
#include <iomanip>
#include <vector>
#include "nu_loss.h"


class Objective : public ceres::CostFunction {
    public:
        Objective(){
            set_num_residuals(1);
            *mutable_parameter_block_sizes() = std::vector<int32_t>(1,2); 
        }
        bool Evaluate(double const* const* parameters,
                              double* residuals,
                              double** jacobians) const {
            residuals[0] = parameters[0][0]*parameters[0][0] - parameters[0][0]*parameters[0][1] + parameters[0][1]*parameters[0][1];
            if (jacobians != nullptr && jacobians[0])
            {
                jacobians[0][0] = 2*parameters[0][0] - parameters[0][1];
                jacobians[0][1] = - parameters[0][0] + 2*parameters[0][1];
            }
            return true;
        }
};

class Constraint : public ceres::CostFunction
{
    public:
        Constraint(){
            set_num_residuals(1);
            *mutable_parameter_block_sizes() = std::vector<int32_t>(1,2);
        }
        bool Evaluate(double const* const* parameters,
                              double* residuals,
                              double** jacobians) const {
            auto g = -parameters[0][1] + 0.75;
            if (g>0)
            {
                residuals[0] = g;
                if (jacobians != nullptr && jacobians[0])
                {
                    jacobians[0][0] = 0;
                    jacobians[0][1] = -1;
                }
                return true;
            }
            residuals[0] = 0.0;
            if (jacobians != nullptr && jacobians[0])
            {
                jacobians[0][0] = 0;
                jacobians[0][1] = 0;
            }
            return true;
        }
};




// MAIN program
int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);
    
    // create decision variables
    std::vector<std::vector<double>> x;
    x.push_back(std::vector<double>{1.8, 1.8});
    
    // Define problem
    ceres::Problem problem;
    problem.AddParameterBlock(x[0].data(), 2);

    //f(x) = x**2 -x*y +y**2
    ceres::CostFunction* cost_function = new Objective;
    problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(nullptr, 2.0, ceres::DO_NOT_TAKE_OWNERSHIP), x[0].data());
    
    // g(x): -y+0.75 <= 0
    ceres::CostFunction* constraint_function = new Constraint;
    ceres::LossFunction* nu_loss = new NuLoss(nullptr, 2*1.0, ceres::DO_NOT_TAKE_OWNERSHIP);
    problem.AddResidualBlock(constraint_function, nu_loss, x[0].data());


    // Set solver parameters
    double param_tolerance = 1e-4;
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    options.max_num_iterations = 100;
    options.parameter_tolerance = 10*param_tolerance;
    ceres::Solver::Summary summary;
    
    // Record x[t] and x[t-1] so for stopping criteria
    Eigen::ArrayXd x_old(2);
    Eigen::ArrayXd x_new(2);
    for (uint i = 0; i<2; ++i)
    {
        memcpy(x_old.data(), x[0].data(), 2*sizeof(double));
    }

    //Inner loop
    while(true)
    {
        // Unconstrained optimization using Ceres
        ceres::Solve(options, &problem, &summary);
        // Print Update or something
        // std::cout << "x: [-1.8, 1.8] -> [" << x[0][0] << ", " << x[0][1] << "]\t nu: " << NuLoss::nu << "\n\n";

        // Copy new data so I can compare it
        for (uint i = 0; i<2; ++i)
        {
            memcpy(x_new.data(), x[0].data(), 2*sizeof(double));
        }

        // Stopping criteria 
        auto temp = x_old - x_new;
        if (temp.abs().any() <= param_tolerance)
        {
            std::cout << "x: [-1.8, 1.8] -> [" << x[0][0] << ", " << x[0][1] << "]\t nu: " << NuLoss::nu << "\n------\n";
            break;
        }
        x_old = x_new;

        // update nu parameter
        NuLoss::nu *=10;
    }
    return 0;
}

