#include "nu_loss.h"

double NuLoss::nu = 1;

NuLoss::NuLoss(const LossFunction* rho, double a, ceres::Ownership ownership)
    : rho_(rho), a_(a), ownership_(ownership)
{
    NuLoss::nu = a;
}


    
void 
NuLoss::Evaluate(double s, double rho[3]) const
{
    if (rho_.get() == nullptr) {
        rho[0] = NuLoss::nu * s;
        rho[1] = NuLoss::nu;
        rho[2] = 0.0;
    }
}