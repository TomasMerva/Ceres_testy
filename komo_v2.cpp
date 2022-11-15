#include "ceres/ceres.h"
#include "glog/logging.h"
#include <iomanip>
#include <vector>

struct KOMO
{
    template <typename T>
    bool operator()(const T* const x_t, const T* const x_t1, const T* const x_t2, T* residual) const {
        // auto feature = (x_t[0] - 2.0*x_t1[0] + x_t2[0]) * (x_t[0] - 2.0*x_t1[0] + x_t2[0]) +
        //                (x_t[1] - 2.0*x_t1[1] + x_t2[1]) * (x_t[1] - 2.0*x_t1[1] + x_t2[1]) +
        //                (x_t[2] - 2.0*x_t1[2] + x_t2[2]) * (x_t[2] - 2.0*x_t1[2] + x_t2[2]) +
        //                (x_t[3] - 2.0*x_t1[3] + x_t2[3]) * (x_t[3] - 2.0*x_t1[3] + x_t2[3]) +
        //                (x_t[4] - 2.0*x_t1[4] + x_t2[4]) * (x_t[4] - 2.0*x_t1[4] + x_t2[4]) +
        //                (x_t[5] - 2.0*x_t1[5] + x_t2[5]) * (x_t[5] - 2.0*x_t1[5] + x_t2[5]) +
        //                (x_t[6] - 2.0*x_t1[6] + x_t2[6]) * (x_t[6] - 2.0*x_t1[6] + x_t2[6]);
        // residual[0] = feature;

        residual[0] = x_t[0] - 2.0*x_t1[0] + x_t2[0];
        residual[1] = x_t[1] - 2.0*x_t1[1] + x_t2[1];
        residual[2] = x_t[2] - 2.0*x_t1[2] + x_t2[2];
        residual[3] = x_t[3] - 2.0*x_t1[3] + x_t2[3];
        residual[4] = x_t[4] - 2.0*x_t1[4] + x_t2[4];
        residual[5] = x_t[5] - 2.0*x_t1[5] + x_t2[5];
        residual[6] = x_t[6] - 2.0*x_t1[6] + x_t2[6];
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
        problem.AddResidualBlock(new ceres::AutoDiffCostFunction<KOMO, 7, 7, 7, 7>(new KOMO), nullptr, x[t], x[t+1], x[t+2]);
    }


    // uint counter=1;
    // while (true)
    // {
    //     // Run the solver!
    //     ceres::Solver::Options options;
    //     options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    //     options.minimizer_progress_to_stdout = true;
    //     options.max_num_iterations = 1;
    //     options.parameter_tolerance = 1e-4;
    //     ceres::Solver::Summary summary;
    //     ceres::Solve(options, &problem, &summary);
    //     std::cout << summary.FullReport() << "\n";

    //     std::cout << "Iteration: " << counter << "\n";
    //     // for (auto xi : x)
    //     // {
    //     //     for (uint i=0; i<num_timestep_var; ++i)
    //     //     {
    //     //         std::cout << xi[i] << "    ";
    //     //     }
    //     //     std::cout << "\n";
    //     // }
    //     counter++;
    //     do 
    //     {
    //     } while (std::cin.get() != '\n');
    // }

    

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
        std::cout << "New timestep" << "\n";
        for (uint i=0; i<num_timestep_var; ++i)
        {
            std::cout << xi[i] << "    ";
        }
        std::cout << "\n";
    }

    // double cost = 0.0;
    // problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

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
