#pragma once
#include "ceres/ceres.h"
#include <Eigen/Dense>

// class Feature_CostFunction : public ceres::CostFunction
// {
//     public:
//         Feature_CostFunction();

//         virtual bool Evaluate(double const* const* parameters,
//                               double* residuals,
//                               double** jacobians) const;

//         void Test_CostFunction();
        
//         double x_vec[7] = {1, 2, 3, 4, 5, 6, 7};

//     private:
//         std::shared_ptr<ceres::Problem> _ceresProblem;

//         const uint _num_joints = 1;
//         const uint _num_timesteps = 7;
// };