#pragma once

#include "ceres/ceres.h"


class NuLoss : public ceres::LossFunction {
    public:
        NuLoss(const LossFunction* rho, double a, ceres::Ownership ownership);
        NuLoss(const NuLoss&) = delete;
        void operator=(const NuLoss&) = delete;
        ~NuLoss() override {
            if (ownership_ == ceres::DO_NOT_TAKE_OWNERSHIP) {
            rho_.release();
            }
        }

        void Evaluate(double s, double rho[3]) const;
        
        static double nu; 

    private:
        std::unique_ptr<const ceres::LossFunction> rho_;
        const ceres::Ownership ownership_;  
        const double a_;
};