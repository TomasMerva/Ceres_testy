#include "ceres/ceres.h"
#include "glog/logging.h"

#include <vector>


class Feature_CostFunction : public ceres::CostFunction
{
    public:
        Feature_CostFunction(){
            set_num_residuals(1);
            *mutable_parameter_block_sizes() = std::vector<int32_t>{2, 2, 2};
        };

        virtual bool Evaluate(double const* const* parameters,
                              double* residuals,
                              double** jacobians) const
        {
            std::cout << parameters[0][0] << " " << parameters[0][1] << " " << parameters[0][2] << "\n";
            std::cout << parameters[1][0] << " " << parameters[1][1] << " " << parameters[1][2] << "\n";
            std::cout << parameters[2][0] << " " << parameters[2][1] << " " << parameters[2][2] << "\n";
            double x = parameters[0][0];
            residuals[0] = 10 - x;
            //     // Compute jacobians
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                jacobians[0][0] = -1;
            }
            return true;
        };

        template <typename T>
        bool operator()(const T* const x, T* residual) const {
            residual[0] = 10.0 - x[0];
            return true;
        }
        

    // private:
    //     std::shared_ptr<ceres::Problem> _ceresProblem;

    //     const uint _num_joints = 1;
    //     const uint _num_timesteps = 7;
};




struct KOMO
{
    template <typename T>
    bool operator()(const T* const x_t, const T* const x_t1, const T* const x_t2, T* residual) const {
        residual[0] = x_t[0] + 10.0 * x_t2[0];
        return true;
    }
};


int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);

    const uint num_timesteps = 4;
    const uint num_x_t = 2;

    double x_t[2] = {1,2};
    double x_t1[2] = {3,4};
    double x_t2[2] = {5,6};
    double x_t3[2] = {7,8};


    ceres::Problem problem;
    problem.AddParameterBlock(x_t, num_x_t);
    problem.AddParameterBlock(x_t1, num_x_t);
    problem.AddParameterBlock(x_t2, num_x_t);
    problem.AddParameterBlock(x_t3, num_x_t);



    problem.AddResidualBlock(new ceres::AutoDiffCostFunction<KOMO, 1, 2, 2, 2>(new KOMO), nullptr, x_t, x_t1, x_t2);
    problem.AddResidualBlock(new ceres::AutoDiffCostFunction<KOMO, 1, 2, 2, 2>(new KOMO), nullptr, x_t1, x_t2, x_t3);

    double cost = 0.0;
    problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    std::cout << "evaluated cost = " << cost << std::endl;


    // Eigen::MatrixXd x_matrix(num_x_t, num_timesteps);
    // // this part is not important
    // x_matrix << 1, 3, 5, 7, 2, 4, 6, 8;
    // std::cout << x_matrix << "\n";


    // std::vector<double *> x;
    // for (uint col=0; col<num_timesteps; ++col)
    // {
    //     x.push_back( new double[num_x_t]{x_matrix(0, col), x_matrix(1, col)} );
    // }

    // ceres::Problem problem;
    // for (auto x_i : x)
    // {
    //     problem.AddParameterBlock(x_i, num_x_t);
    // }

    // // upload all blocks to cost function
    // ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<KOMO, 1, num_timesteps>(new KOMO);
    // problem.AddResidualBlock(cost_function, nullptr, x); 

    // double cost = 0.0;
    // problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

//     std::vector<double> x_t{1,2};
//     std::vector<double> x_t1{3,4};
//     std::vector<double> x_t2{5,4};


    
//     // double x_t[2] = {1, 2};
//     // double x_t1[2] = {3, 4};
//     // double x_t2[2] = {5, 6};

//     // std::vector<double> x{1, 2, 3, 4};
//     // double x[3] = {1, 2, 5};

//     // std::vector<double*> x{x_t, x_t1, x_t2};

//     ceres::Problem problem;
//     // problem.AddParameterBlock(x, 3); // not necessary
//     problem.AddParameterBlock(x_t.data(), 2);
//     problem.AddParameterBlock(x_t1.data(), 2);
//     problem.AddParameterBlock(x_t2.data(), 2);
//     // std::vector<double *> x{x_t, x_t1, x_t2};
//     std::vector<std::vector<double> *> x{x_t.data(), x_t1.data(), x_t2.data()};


//     // 

//     // ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<KOMO_Functor, 1, 3>(new KOMO_Functor);
//     ceres::CostFunction* cost_function = new Feature_CostFunction;

// //   std::vector<int32_t>* mutable_parameter_block_sizes() {
// //     return &parameter_block_sizes_;
// //   }


//     problem.AddResidualBlock(cost_function, nullptr, x);  //num parameters block is defined by x, since x is a vector of 3 pointers, then it expects 3 blocks

//     double cost = 0.0;
//     problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
//     std::cout << "evaluated cost = " << cost << std::endl;

    // // problem.EvaluateResidualBlock(id, true, &cost,  nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

    // for (uint t=0; t<(num_timesteps-2); ++t)
    // {
    //     std::cout << t << "\n";
    //     problem.AddParameterBlock(&x[t], 3);
    // }



    // problem.AddParameterBlock(x.data(), x.size());
    // auto id = problem.AddResidualBlock(new ceres::AutoDiffCostFunction<t1, 1, 4>(new t1), nullptr, x.data());

    // problem.AddResidualBlock(
    //   new ceres::AutoDiffCostFunction<F1, 1, 1, 1>(new F1), nullptr, &x1, &x2);

    // problem.AddResidualBlock(
    //   new ceres::AutoDiffCostFunction<F2, 1, 1, 1>(new F2), nullptr, &x3, &x4);

    
    // problem.SetParameterLowerBound(x.data(), 0, -1.2116);
    // problem.SetParameterUpperBound(x.data(), 0, -1.2116);
    // problem.SetParameterLowerBound(x.data(), 3, -0.507423);
    // problem.SetParameterUpperBound(x.data(), 3, -0.507423);


    // double cost = 0.0;
    // problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, nullptr, nullptr, nullptr);
    
    // // problem.EvaluateResidualBlock(id, true, &cost,  nullptr, nullptr);
    // std::cout << "evaluated cost = " << cost << std::endl;

    // Run the solver!
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << "\n";
    // std::cout << "x : " << initial_x << " -> " << x[0] << "\n";

    return 0;

}