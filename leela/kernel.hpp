//
//  kernel.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef kernel_hpp
#define kernel_hpp

#include <iostream>
#include <cmath>
#include "data.h"

double kernel_wendland_alpha_1d();
double kernel_wendland_alpha_2d();
double kernel_wendland_alpha_3d();
double kernel_wendland(const vec3<double> &r, double alpha, double h, int dim);
double kernel_wendland_dash(const vec3<double> &r, double alpha, double h, int dim);
vec3<double> kernel_e(const vec3<double> &r);

#endif /* kernel_hpp */
