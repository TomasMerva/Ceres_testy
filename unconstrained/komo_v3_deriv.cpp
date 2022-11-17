#include "ceres/ceres.h"
#include "glog/logging.h"
#include <iomanip>
#include <vector>

class KOMO_FeatureCost : public ceres::CostFunction {
    public:
        KOMO_FeatureCost() {
            set_num_residuals(7);
            *mutable_parameter_block_sizes() = std::vector<int32_t>(3, 7);
            // *mutable_parameter_block_sizes() = std::vector<int32_t>{7,7,7};
        }
        bool Evaluate(double const* const* parameters,
                              double* residuals,
                              double** jacobians) const {
            // residuals[0] = 10 - parameters[0][0];
            residuals[0] = parameters[0][0] - 2.0*parameters[1][0] + parameters[2][0];
            residuals[1] = parameters[0][1] - 2.0*parameters[1][1] + parameters[2][1];
            residuals[2] = parameters[0][2] - 2.0*parameters[1][2] + parameters[2][2];
            residuals[3] = parameters[0][3] - 2.0*parameters[1][3] + parameters[2][3];
            residuals[4] = parameters[0][4] - 2.0*parameters[1][4] + parameters[2][4];
            residuals[5] = parameters[0][5] - 2.0*parameters[1][5] + parameters[2][5];
            residuals[6] = parameters[0][6] - 2.0*parameters[1][6] + parameters[2][6];

            // Compute the Jacobian if asked for.
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                // jacobians[id_parameterblocku][id_res*block_size + id_var]
                memset(jacobians[0], 0, 7*7*sizeof(jacobians[0]));
                memset(jacobians[1], 0, 7*7*sizeof(jacobians[1]));
                memset(jacobians[2], 0, 7*7*sizeof(jacobians[2]));
                // res[0]
                jacobians[0][0] = 1.0;
                jacobians[1][0] = -2.0;
                jacobians[2][0] = 1.0;
                // res[1]
                jacobians[0][8] = 1.0;
                jacobians[1][8] = -2.0;
                jacobians[2][8] = 1.0;
                // res[2]
                jacobians[0][16] = 1.0;
                jacobians[1][16] = -2.0;
                jacobians[2][16] = 1.0;
                // res[3]
                jacobians[0][24] = 1.0;
                jacobians[1][24] = -2.0;
                jacobians[2][24] = 1.0;
                // res[4]
                jacobians[0][32] = 1.0;
                jacobians[1][32] = -2.0;
                jacobians[2][32] = 1.0;
                // res[5]
                jacobians[0][40] = 1.0;
                jacobians[1][40] = -2.0;
                jacobians[2][40] = 1.0;
                // res[6]
                jacobians[0][48] = 1.0;
                jacobians[1][48] = -2.0;
                jacobians[2][48] = 1.0;
            }
            return true;
        }
};

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in);

int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);
    uint num_timesteps = 20;
    uint num_timestep_var = 7;

    
    // generate q linspace
    std::vector<double> q_start{-1.2116, 0.356153, 0.244033, -2.28525, -1.68591, 2.35773, 4.035e-08};
    std::vector<double> q_end{-1.06086, 0.436372, 0.461309, -2.32378, -1.67167, 2.46085, 4.035e-08};
    auto q0 = linspace(q_start[0], q_end[0], num_timesteps);
    auto q1 = linspace(q_start[1], q_end[1], num_timesteps);
    auto q2 = linspace(q_start[2], q_end[2], num_timesteps);
    auto q3 = linspace(q_start[3], q_end[3], num_timesteps);
    auto q4 = linspace(q_start[4], q_end[4], num_timesteps);
    auto q5 = linspace(q_start[5], q_end[5], num_timesteps);
    auto q6 = linspace(q_start[6], q_end[6], num_timesteps);


    // create vector of pointers to arrays representing timesteps
    std::vector<double *> x;
    for (uint timestep=0; timestep<num_timesteps; ++timestep)
    {
        // x.push_back( new double[num_timestep_var]{q0[timestep], q1[timestep], q2[timestep], q3[timestep], q4[timestep], q5[timestep], q6[timestep]} );
        x.push_back( new double[num_timestep_var]{q0[0], q1[0], q2[0], q3[0], q4[0], q5[0], q6[0]} );
    }
    // end state
    x[19][0] = q_end[0];
    x[19][1] = q_end[1];
    x[19][2] = q_end[2];
    x[19][3] = q_end[3];
    x[19][4] = q_end[4];
    x[19][5] = q_end[5];
    x[19][6] = q_end[6];


    ceres::Problem problem;
    for (uint t=0; t<num_timesteps; ++t)
    {
        problem.AddParameterBlock(x[t], num_timestep_var);
        // q1
        problem.SetParameterLowerBound(x[t], 0, -2.8973);
        problem.SetParameterUpperBound(x[t], 0, 2.8973);
        // q2
        problem.SetParameterLowerBound(x[t], 1, -1.7628);
        problem.SetParameterUpperBound(x[t], 1, 1.7628);
        // q3
        problem.SetParameterLowerBound(x[t], 2, -2.8973);
        problem.SetParameterUpperBound(x[t], 2, 2.8973);
        // q4
        problem.SetParameterLowerBound(x[t], 0, -3.0718);
        problem.SetParameterUpperBound(x[t], 0, -0.0698);
        // q5
        problem.SetParameterLowerBound(x[t], 1, -2.8973);
        problem.SetParameterUpperBound(x[t], 1, 2.8973);
        // q6
        problem.SetParameterLowerBound(x[t], 2, -0.0175);
        problem.SetParameterUpperBound(x[t], 2, 3.7525);
        // q7
        problem.SetParameterLowerBound(x[t], 2, -2.8973);
        problem.SetParameterUpperBound(x[t], 2, 2.8973);
    }
    for (uint i=0; i<num_timestep_var; ++i)
    {
        // Init state boundaries
        problem.SetParameterLowerBound(x[0], i, q_start[i]-0.0001);
        problem.SetParameterUpperBound(x[0], i, q_start[i]+0.0001);
         // End state boundaries
        problem.SetParameterLowerBound(x[num_timesteps-1], i, q_end[i]-0.0001);
        problem.SetParameterUpperBound(x[num_timesteps-1], i, q_end[i]+0.0001);
    }

    for (uint t=0; t<(num_timesteps-2); ++t)
    {
        ceres::CostFunction* cost_function = new KOMO_FeatureCost;
        problem.AddResidualBlock(cost_function, nullptr, x[t], x[t+1], x[t+2]);
    }

    // Run the solver!
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 1000;
    options.parameter_tolerance = 1e-4;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";

    for (auto xi : x)
    {
        for (uint i=0; i<num_timestep_var; ++i)
        {
            std::cout << xi[i] << "    ";
        }
        std::cout << "\n";
    }
    return 0;
}

///////////////////////////////////////////////////////

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}
