//
//  solve.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef solve_hpp
#define solve_hpp

#include <iostream>
#include <cmath>
#include "data.h"
#include "utils.hpp"

void solve(data &d);
void solve_timestep(data &d);
double solve_timestep_cfl(data &d);
double solve_timestep_viscosity(data &d);
double solve_timestep_surfacetension(data &d);
void solve_timestep_output(data &d);
void solve_continuity(data &d);
void solve_continuity_density(data &d);
void solve_continuity_smoothing(data &d);
void solve_eos(data &d);
void solve_eos_pressure(data &d);
void solve_momentum(data &d);
void solve_momentum_acceleration(data &d);
std::vector<vec3<double>> solve_momentum_acceleration_pressure(data &d);
std::vector<vec3<double>> solve_momentum_acceleration_viscosity(data &d);
std::vector<vec3<double>> solve_momentum_acceleration_surfacetension(data &d);
std::vector<vec3<double>> solve_momentum_acceleration_gravity(data &d);
std::vector<vec3<double>> solve_momentum_acceleration_tensileinstability(data &d);
void solve_verlet(data &d);
void solve_verlet_velocity(data &d);
void solve_verlet_position(data &d);

#endif /* solve_hpp */
