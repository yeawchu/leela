//
//  kernel.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "kernel.hpp"

double kernel_wendland_alpha_1d()
{
    return 3 / 4;
}

double kernel_wendland_alpha_2d()
{
    return 7 / (4 * M_PI);
}

double kernel_wendland_alpha_3d()
{
    return 21 / (16 * M_PI);
}

double kernel_wendland(const vec3<double> &r, double alpha, double h, int dim)
{
    // wendland kernel, h_kappa = 2
    // for 0 <= s <= 2, W = (W_alpha / h^n) (1 - s / 2)^4 (1 + 2 s)
    // for 2 <  s     , W = 0
    double s = mag(r) / h;
    if (s <= 2.0) return alpha / std::pow(h, static_cast<double>(dim)) * std::pow(1.0 - 0.5 * s, 4.0) * (1.0 + 2.0 * s);
    else return 0;
}

double kernel_wendland_dash(const vec3<double> &r, double alpha, double h, int dim)
{
    // wendland kernel, h_kappa = 2
    // W' = (W_alpha / h^(n + 1)) dW/ds
    // for 0 <= s <= 2, W' = (W_alpha / h^(n + 1)) (-5 s(1 - s / 2)^3)
    // for 2 < s      , W' = 0
    double s = mag(r) / h;
    if (s <= 2.0) return  alpha / std::pow(h, static_cast<double>(dim) + 1.0) * (-5.0 * s * std::pow(1.0 - 0.5 * s, 3.0));
    else return 0;
}

vec3<double> kernel_e(const vec3<double> &r)
{
    // unit vector r
    // r/|r|
    if (mag(r) != 0) return normalise(r);
    else return vec3<double> (0, 0, 0);
}