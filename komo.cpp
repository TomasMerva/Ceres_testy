#include "komo.h"


// Feature_CostFunction::Feature_CostFunction()
// {
//     _ceresProblem = std::make_shared<ceres::Problem>();

//     for (uint t=0; t<_num_timesteps; ++t)
//     {
//         // _ceresProblem->AddParameterBlock(x_vec)
//     }
// }

// bool 
// Feature_CostFunction::Evaluate
//     (const double* const* parameters, double* residuals, double** jacobians) const 
// {
//     // const double* p -> points to constant but it can change its address
//     // double* const p -> pointer cannot change its address but it does not point to const variable
//     // const double* const p -> const pointer to const variable

//     // Compute jacobians
//     if (jacobians != nullptr && jacobians[0] != nullptr) {
//       jacobians[0][0] = -1;
//     }
//     return true;
// }

// void
// Feature_CostFunction::Test_CostFunction()
// {
//  // Evaluate();
// }

