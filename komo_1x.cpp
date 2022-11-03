#include "ceres/ceres.h"
#include "glog/logging.h"

#include <vector>

// struct Objective
// {
//     template <typename T>
//     bool operator()(const T* parameters, T* residual) const 
//     {
//         const T x = parameters[0];
//         for (uint i=0; i<7; ++i)
//         {
//             const T temp = parameters[i];
//             std::cout << temp << "    ";
//         }
//         std::cout << "\n\n\n";
//         residual[0] = (7.0-x);
//         return true;
//     }
// };

struct KOMO
{
    template <typename T>
    bool operator()(const T* parameters, T* residual) const 
    {
        const T x = parameters[0];
        residual[0] = (7.0-x);
        return true;
    }
};



int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);

    double initial_x = 1.0;

    std::vector<double> x{-1.2116, -1.2116, -1.2116, -0.507423};

    ceres::Problem problem;
    problem.AddParameterBlock(x.data(), 3);
    auto id = problem.AddResidualBlock(new ceres::AutoDiffCostFunction<KOMO, 1, 1>(new KOMO), nullptr, x.data());
    problem.SetParameterLowerBound(x.data(), 0, -1.2116);
    problem.SetParameterUpperBound(x.data(), 0, -1.2116);
    problem.SetParameterLowerBound(x.data(), 3, -0.507423);
    problem.SetParameterUpperBound(x.data(), 3, -0.507423);



    // double cost = 0.0;
    // problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    
    // problem.EvaluateResidualBlock(id, true, &cost,  nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

    // Run the solver!
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << "\n";
    std::cout << "x : " << initial_x << " -> " << x[0] << "\n";

    return 0;

}