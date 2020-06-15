//
//  solve.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "solve.hpp"

// main solver
void solve(data &d)
{
    solve_timestep(d);
    solve_continuity(d);
    solve_eos(d);
    solve_momentum(d);
    solve_verlet(d);
}

// determine time step size
void solve_timestep(data &d)
{
    double dt_cfl = solve_timestep_cfl(d);
    double dt_vis = solve_timestep_viscosity(d);
    //double dt_st = solve_timestep_surfacetension(d);
    
    d.time.dt = vmin(dt_cfl, dt_vis);         // find the minimum/smallest time step to satisfy all particles
    //d.time.dt = vmin(dt_cfl, dt_vis, dt_st);         // find the minimum/smallest time step to satisfy all particles

    solve_timestep_output(d);
}

// cfl criteria
double solve_timestep_cfl(data &d)
{
    double dt_cfl_min = d.time.dt_max;
    double dt_cfl = d.time.dt_c_cfl * d.constant.h / d.constant.c_0;
    dt_cfl_min = std::min(dt_cfl_min, dt_cfl);
    return dt_cfl_min;
}

// viscous criteria
double solve_timestep_viscosity(data &d)
{
    double dt_vis_min = d.time.dt_max;
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        double dt_vis = d.time.dt_c_viscosity * d.constant.h * d.constant.h * d.fluid.rho[i] / d.fluid.mu[i];
        dt_vis_min = std::min(dt_vis_min, dt_vis);
    }
    return dt_vis_min;
}

double solve_timestep_surfacetension(data &d)
{
    double dt_st_min = d.time.dt_max;
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        double dt_st = d.time.dt_c_surfacetension * std::sqrt(d.constant.rho_0 * d.constant.h * d.constant.h * d.constant.h / d.fluid.gamma[i]);
        dt_st_min = std::min(dt_st_min, dt_st);
    }
    return dt_st_min;
}

// Adjust timestep to match output interval
void solve_timestep_output(data &d)
{
    if ((d.time.t + d.time.dt) > static_cast<double>(d.io.itr + 1) * d.time.dt_output) // resize dt if exceed output_itr
    {
        ++d.io.itr;                                     // increase output iterator
        d.time.dt = static_cast<double>(d.io.itr) * d.time.dt_output - d.time.t;
    }
    else if ((d.time.t + 1.1 * d.time.dt) > static_cast<double>(d.io.itr + 1) * d.time.dt_output) // introduced a factor of 1.1 to capture float precision problems when t+dt is marginally > (itr+1)*output_itr.  If not, the very small computed dt will be potential for errors as anything divide by very small dt -> 0 is going to result in very large values!!!
    {
        d.time.dt = 0.5 * (static_cast<double>(d.io.itr + 1) * d.time.dt_output - d.time.t);    // used a factor or 0.5 here to average out the time step.
    }
}

void solve_continuity(data &d)
{
    solve_continuity_density(d);
    solve_continuity_smoothing(d);
}

// continuum density
void solve_continuity_density(data &d)
{
    std::vector<double> sum(d.fluid.pos.size());
    
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        sum[i] += d.fluid.m[j] * (dot((d.fluid.vel[i] - d.fluid.vel[j]), d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a]));
    }
    
    for (long a = 0, a_size = d.nlist.fb.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.fb.i[a];
        long j = d.nlist.fb.j[a];
        vec3<double> vel_k = 2.0 * d.boundary.vel[j] - d.fluid.vel[i];
        sum[i] += d.fluid.m[i] * (dot((d.fluid.vel[i] - vel_k), d.nlist.fb.w_e[a] * d.nlist.fb.w_dash[a]));
    }

    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        d.fluid.rho_n[i] = d.fluid.rho_o[i] + sum[i] * (d.time.dt + d.time.dt_o);
    }
}

// density gravity smoothing - similar to xsph - see Violeau pg 488
void solve_continuity_smoothing(data &d)
{
    std::vector<double> sum(d.fluid.pos.size());
    
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        sum[i] += d.fluid.m[j] * (d.fluid.rho[i] - d.fluid.rho[j]) / (0.5 * (d.fluid.rho[i] + d.fluid.rho[j])) * d.nlist.ff.w[a];
    }
    // Note that assessing boundary particles give 0 because rho[i] - rho[j] = 0
    double epsilon_rho = 0.05;      // xpsh-like factor for density
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        d.fluid.rho[i] = d.fluid.rho[i] - epsilon_rho * sum[i];
    }
}

// pressure
void solve_eos(data &d)
{
    solve_eos_pressure(d);
}

// pressure using tait's equation
void solve_eos_pressure(data &d)
{
    double gamma = 7.0;
    double b = d.constant.rho_0 * d.constant.c_0 * d.constant.c_0 / gamma;
    
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        // tait's equation
        // variation of density remains small and volume of fluid is generally well conserved
        // limitations: isothermal conditions
        d.fluid.p[i] = b * (std::pow((d.fluid.rho[i] / d.constant.rho_0), gamma) - 1.0);
    }
}

void solve_momentum(data &d)
{
    solve_momentum_acceleration(d);
}

// acceleration
void solve_momentum_acceleration(data &d)
{
    auto acc_p = solve_momentum_acceleration_pressure(d);
    auto acc_mu = solve_momentum_acceleration_viscosity(d);
    //auto acc_st = solve_momentum_acceleration_surfacetension(d);
    auto acc_g = solve_momentum_acceleration_gravity(d);

    //auto acc_ti = solve_momentum_acceleration_tensileinstability(d);

    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        d.fluid.acc[i] = acc_p[i] + acc_mu[i] + acc_g[i];                  // total acceleration
        //d.fluid.acc[i] = acc_p[i] + acc_mu[i] + acc_st[i] + acc_g[i];                  // total acceleration
        //d.fluid.acc[i] = acc_p[i] + acc_mu[i] + acc_st[i] + acc_g[i] + acc_ti[i];                  // total acceleration
    }
}

std::vector<vec3<double>> solve_momentum_acceleration_pressure(data &d)
{
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        // pressure
        sum[i] += - d.fluid.m[j] * (d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i]) + d.fluid.p[j] / (d.fluid.rho[j] * d.fluid.rho[j])) * d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a];
    }
    
    for (long a = 0, a_size = d.nlist.fb.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.fb.i[a];
        
        // pressure
        sum[i] += - d.fluid.m[i] * (d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i]) + d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i])) * d.nlist.fb.w_e[a] * d.nlist.fb.w_dash[a];
    }
    return sum;
}

std::vector<vec3<double>> solve_momentum_acceleration_viscosity(data &d)
{
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        // viscosity
        double mag_r_ij = mag(d.fluid.pos[i] - d.fluid.pos[j]);
        double eta = 0.1 * d.constant.h;
        vec3<double> vel_ij = d.fluid.vel[i] - d.fluid.vel[j];
        sum[i] += d.fluid.m[j] * (d.fluid.mu[i] + d.fluid.mu[j]) / (2.0 * d.fluid.rho[i] * d.fluid.rho[j]) * mag_r_ij / (mag_r_ij * mag_r_ij + eta * eta) * (static_cast<double>(d.constant.dim + 2) * dot(vel_ij, d.nlist.ff.w_e[a]) * d.nlist.ff.w_e[a] + vel_ij) * d.nlist.ff.w_dash[a];
    }
    
    for (long a = 0, a_size = d.nlist.fb.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.fb.i[a];
        long j = d.nlist.fb.j[a];
        vec3<double> vel_k = 2.0 * d.boundary.vel[j] - d.fluid.vel[i];
        vec3<double> pos_k = 2.0 * d.boundary.pos[j] - d.fluid.pos[i];
        
        // viscosity
        double mag_r_ij = mag(d.fluid.pos[i] - pos_k);
        double eta = 0.1 * d.constant.h;
        vec3<double> vel_ij = d.fluid.vel[i] - vel_k;
        sum[i] += d.fluid.m[i] * (d.fluid.mu[i] + d.fluid.mu[i]) / (2.0 * d.fluid.rho[i] * d.fluid.rho[i]) * mag_r_ij / (mag_r_ij * mag_r_ij + eta * eta) * (static_cast<double>(d.constant.dim + 2) * dot(vel_ij, d.nlist.fb.w_e[a]) * d.nlist.fb.w_e[a] + vel_ij) * d.nlist.fb.w_dash[a];
    }

    return sum;
}

std::vector<vec3<double>> solve_momentum_acceleration_surfacetension(data &d)
{
   /*
     std::vector<double> c(d.fluid.pos.size());
     
     // calculate colour function
     for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
     {
     long i = d.nlist.ff.i[a];
     long j = d.nlist.ff.j[a];
     
     c[i] += d.fluid.m[j] / d.fluid.rho[j] * d.nlist.ff.w[a];
     }

    std::vector<vec3<double>> n(d.fluid.pos.size());
    std::vector<vec3<double>> n_unit(d.fluid.pos.size());
    
    // calculate normal at interface
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        //n[i] += d.fluid.m[j] / d.fluid.rho[j] * d.nlist.ff.w_dash[a] * d.nlist.ff.w_e[a];
        n[i] += d.fluid.rho[i] * d.fluid.m[j] * (c[i] / (d.fluid.rho[i] * d.fluid.rho[i]) + c[j] / (d.fluid.rho[j] * d.fluid.rho[j])) * d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a];
    }
    
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        n_unit[i] = normalise(n[i]);
    }
    
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    // surface tension
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        auto n_ij = n_unit[i] - n_unit[j];
                
        sum[i] += d.fluid.gamma[i] / d.fluid.rho[i] * d.fluid.m[j] / d.fluid.rho[j] * dot(n_ij, d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a]) * n[j];
    }
*/

    std::vector<double> c(d.fluid.pos.size());
    
    // calculate colour function
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        c[i] += d.fluid.m[j] / d.fluid.rho[j] * d.nlist.ff.w[a];
    }
    
    std::vector<vec3<double>> n(d.fluid.pos.size());
    std::vector<vec3<double>> n_unit(d.fluid.pos.size());
    
    // calculate normal at interface
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        //n[i] += d.fluid.m[j] / d.fluid.rho[j] * d.nlist.ff.w_dash[a] * d.nlist.ff.w_e[a];
        n[i] += d.fluid.rho[i] * d.fluid.m[j] * (c[i] / (d.fluid.rho[i] * d.fluid.rho[i]) + c[j] / (d.fluid.rho[j] * d.fluid.rho[j])) * d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a];
    }
    
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        n_unit[i] = normalise(n[i]);
    }
    
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    // surface tension
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        auto n_ij = n_unit[i] - n_unit[j];

        sum[i] += d.fluid.gamma[i] / d.fluid.rho[i] * d.fluid.m[j] / d.fluid.rho[j] * dot(n_ij, d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a]) * n[j];
    }

    return sum;
}

std::vector<vec3<double>> solve_momentum_acceleration_gravity(data &d)
{
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        sum[i] = d.constant.g;                                       // gravity acceleration
    }
    
    return sum;
}

// tensile instability
std::vector<vec3<double>> solve_momentum_acceleration_tensileinstability(data &d)
{
    std::vector<vec3<double>> sum(d.fluid.pos.size());
    
    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        
        double ti_coef_i;
        double ti_coef_j;
        if (d.fluid.p[i] > 0) ti_coef_i = 0.001; else ti_coef_i = -0.1;
        if (d.fluid.p[j] > 0) ti_coef_j = 0.001; else ti_coef_j = -0.1;
        double r = 0.5 * d.constant.h;
        vec3<double> v_r = {std::sqrt(r * r / 3), std::sqrt(r * r / 3), std::sqrt(r * r / 3)};  // Recommended fix by Guillermo to avoid dimensionality bug other than 3D
        // vec3<double> v_r = {std::sqrt(r * r / d.constant.dim), std::sqrt(r * r / d.constant.dim), std::sqrt(r * r / d.constant.dim)};
        double w_q_min = d.kernel.func_w(v_r, d.kernel.alpha, d.constant.h, d.constant.dim);
        double psi = std::pow(d.nlist.ff.w[a] / w_q_min, 4.0);
        
        sum[i] += - d.fluid.m[j] * (d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i]) * ti_coef_i * psi + d.fluid.p[j] / (d.fluid.rho[j] * d.fluid.rho[j]) * ti_coef_j * psi) * d.nlist.ff.w_e[a] * d.nlist.ff.w_dash[a];
    }
    
    for (long a = 0, a_size = d.nlist.fb.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.fb.i[a];
        
        double ti_coef_i;
        if (d.fluid.p[i] > 0) ti_coef_i = 0.001; else ti_coef_i = -0.1;
        double r = 0.5 * d.constant.h;
        vec3<double> v_r = {std::sqrt(r * r / d.constant.dim), std::sqrt(r * r / d.constant.dim), std::sqrt(r * r / d.constant.dim)};
        double w_q_min = d.kernel.func_w(v_r, d.kernel.alpha, d.constant.h, d.constant.dim);
        double psi = std::pow(d.nlist.fb.w[a] / w_q_min, 4.0);

        sum[i] += - d.fluid.m[i] * (d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i]) * ti_coef_i * psi + d.fluid.p[i] / (d.fluid.rho[i] * d.fluid.rho[i]) * ti_coef_i * psi) * d.nlist.fb.w_e[a] * d.nlist.fb.w_dash[a];
    }
    return sum;
}

// stormer-verlet method
void solve_verlet(data &d)
{
    // note we are working on time n and not n+1
    // acceleration with non-constant time steps:
    // x(n+1) = x(n) + (x(n) - x(n-1))(dt(n)/dt(n-1)) + a(n)0.5(dt(n) + dt(n-1))dt(n)
    // v(n+1) = v(n-1) + a(n)(dt(n-1) + dt(n))
    solve_verlet_position(d);
    solve_verlet_velocity(d);
}

// position
void solve_verlet_position(data &d)
{
    double dt_ratio = d.time.dt / d.time.dt_o;
    double dt_ohalf = 0.5 * (d.time.dt + d.time.dt_o);
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        d.fluid.pos_n[i] = d.fluid.pos[i] + (d.fluid.pos[i] - d.fluid.pos_o[i]) * dt_ratio + d.fluid.acc[i] * dt_ohalf * d.time.dt;
    }
}

// velocity
void solve_verlet_velocity(data &d)
{
    // adapted from heun's method or explicit trapizoidal rule
    // u(n+1) = 2 (x(n+1) - x(n)) / dt(n) - u(n)
    for (long i = 0, i_size = d.fluid.pos.size(); i != i_size; ++i)
    {
        d.fluid.vel_n[i] = 2.0 * (d.fluid.pos_n[i] - d.fluid.pos[i]) / d.time.dt - d.fluid.vel[i];
        //d.fluid.vel_n[i] = 2.0 * (d.fluid.pos_n[i] - d.fluid.pos_o[i]) / (d.time.dt + d.time.dt_o) - d.fluid.vel_o[i];
    }
}
