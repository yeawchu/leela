//
//  data.h
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef data_h
#define data_h

#include <vector>
#include <functional>
#include "vec.h"

struct neighbour_list_parameters
{
    std::vector<vec3<double>> w_e;               // (determine live) kernel gradient
    
    std::vector<double> w;                         // (determine live) kernel
    std::vector<double> w_dash;                         // (determine live) kernel
    
    std::vector<long> i;         // (determine live) particle fluid to fluid neighbour list ids
    std::vector<long> j;         // (determine live) particle fluid to fluid neighbour list ids
};

struct neighbour_list
{
    neighbour_list_parameters ff;
    neighbour_list_parameters fb;
};

struct fluid_particle
{
    std::vector<vec3<double>> pos_n;              // (determine live) particle position at n+1
    std::vector<vec3<double>> pos;                // particle position at n
    std::vector<vec3<double>> pos_o;              // particle position at n-1
    std::vector<vec3<double>> vel_n;              // (determine live) particle velocity at n+1
    std::vector<vec3<double>> vel;                // particle velocity at n
    std::vector<vec3<double>> vel_o;              // particle velocity at n-1
    std::vector<vec3<double>> acc;                // particle acceleration at n

    std::vector<double> p;              // particle pressure at n
    std::vector<double> rho_n;          // (determine live) particle density at n+1
    std::vector<double> rho;            // particle density at n
    std::vector<double> rho_o;          // particle density at n-1
    std::vector<double> m;              // particle mass
    std::vector<double> mu;             // particle viscosity
    std::vector<double> gamma;             // particle surface tension
};

struct boundary_particle
{
    std::vector<vec3<double>> pos;           // boundary wall particle position at n
    std::vector<vec3<double>> vel;           // boundary wall particle velocity at n
};

struct kernel_parameters
{
    double alpha;             // (determine live) kernel alpha for dimensional renormalisation

    std::string type;             // kernel type

    std::function<double(const vec3<double> &, double, double, int)> func_w;                  // (determine live) w
    std::function<double(const vec3<double> &, double, double, int)> func_w_dash;                  // (determine live) w'
    std::function<vec3<double>(const vec3<double> &)> func_w_e;        // (determine live) unit vector
};

struct time_parameters
{
    double t_begin;             // start time
    double t_end;               // end time
    double t;                   // (determine live) current time
    double dt;                  // (determine live) time step at n
    double dt_o;                // (determine live) time step at n-1
    double dt_c_cfl;            // time step cfl coefficient
    double dt_c_viscosity;            // time step viscosity coefficient
    double dt_c_surfacetension;            // time step surface tension coefficient
    double dt_max;              // maximum time step
    double dt_output;           // time step for outputs
};

struct constant_parameters
{
    vec3<double> g;             // gravity
    
    double rho_0;              // reference density
    double c_0;                // reference sound speed

    double h_eta;                 // expansion ratio for smoothing length to capture sufficient neighbouring particles
    double h_kappa;               // support domain to smoothing length ratio
    double h;                   // (determine live) constant smoothing length for each particle

    double delta_s;               // default particle spacing

    int dim;                        // dimensionality of the problem
};

struct io_parameters
{
    std::string core_base;      // datafile basename
    std::string core_ext;       // datafile extension
    std::string vtk_base_vtu;           // vtu basename
    std::string vtk_ext_vtu;            // vtu extension
    std::string vtk_base_pvd;           // pvd basename
    std::string vtk_ext_pvd;            // pvd extension

    long itr;                       // (determined live) output file iterator

    int precision;               // input/output precision
};

struct data
{
    neighbour_list nlist;
    fluid_particle fluid;
    boundary_particle boundary;
    kernel_parameters kernel;
    time_parameters time;
    constant_parameters constant;
    io_parameters io;
};

#endif /* data_h */
