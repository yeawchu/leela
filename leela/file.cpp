//
//  file.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "file.hpp"

// read init file
void file_config_read(data &d, const std::string &filename)
{
    std::ifstream ifs;
    
    ifs.open (filename, std::ifstream::in);
    if (ifs.is_open())                                  // Check if file is open
    {
        long counter = 1;                               // line number counter
        while (ifs.good())
        {
            std::string line;
            while (std::getline(ifs, line))
            {
                file_config_read_parse(d, line, counter);
                ++ counter;
            }
        }
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
}

// read parse lines in init file
void file_config_read_parse(data &d, const std::string &line, long counter)
{
    std::vector<std::string> v = tokenise<std::string>(line, " \t");    // tokenise using space and tab
    for (long i = 0, i_size=v.size(); i != i_size; ++i)
    {                                                                   // check token is
        if ((v[i] == "#") ||                                            // #
            (v[i] == ";") ||                                            // ;
            (v.empty()) ||                                              // empty line
            (v[i][0] == '#') ||                                         // begin with #
            (v[i][0] == ';'))                                           // begin with ;
        {
            v.erase(v.begin() + i, v.begin() + i_size);                 // erase token and everthing after it
            break;                                                      // break from for loop
        }
        std::transform(v[i].begin(), v[i].end(), v[i].begin(), ::tolower); // converts to lower case
    }
    file_config_read_label(d, v, counter);
}

// read options
void file_config_read_label(data &d, const std::vector<std::string> &v, long counter)
{
    for (long i = 0, i_size = v.size(); i < i_size; ++i)
    {
        // constants
        if (v[i] == "constant.dim") { file_read_label_var(v, d.constant.dim, i); continue; }
        if (v[i] == "constant.g") { file_read_label_vec(v, d.constant.g, i); continue; }
        
        // reference parameters
        if (v[i] == "constant.rho_0") { file_read_label_var(v, d.constant.rho_0, i); continue; }
        if (v[i] == "constant.c_0") { file_read_label_var(v, d.constant.c_0, i); continue; }
        
        // sph kernel
        if (v[i] == "constant.h_eta") { file_read_label_var(v, d.constant.h_eta, i); continue; }
        if (v[i] == "constant.h_kappa") { file_read_label_var(v, d.constant.h_kappa, i); continue; }
        if (v[i] == "kernel.type") { file_read_label_var(v, d.kernel.type, i); continue; }
        
        // time
        if (v[i] == "time.t_begin") { file_read_label_var(v, d.time.t_begin, i); continue; }
        if (v[i] == "time.t_end") { file_read_label_var(v, d.time.t_end, i); continue; }
        if (v[i] == "time.dt_max") { file_read_label_var(v, d.time.dt_max, i); continue; }
        if (v[i] == "time.dt_c_cfl") { file_read_label_var(v, d.time.dt_c_cfl, i); continue; }
        if (v[i] == "time.dt_c_viscosity") { file_read_label_var(v, d.time.dt_c_viscosity, i); continue; }
        if (v[i] == "time.dt_c_surfacetension") { file_read_label_var(v, d.time.dt_c_surfacetension, i); continue; }
        if (v[i] == "time.dt_output") { file_read_label_var(v, d.time.dt_output, i); continue; }
        
        // io
        if (v[i] == "io.core_base") { file_read_label_var(v, d.io.core_base, i); continue; }
        if (v[i] == "io.core_ext") { file_read_label_var(v, d.io.core_ext, i); continue; }
        if (v[i] == "io.vtk_base_vtu") { file_read_label_var(v, d.io.vtk_base_vtu, i); continue; }
        if (v[i] == "io.vtk_ext_vtu") { file_read_label_var(v, d.io.vtk_ext_vtu, i); continue; }
        if (v[i] == "io.vtk_base_pvd") { file_read_label_var(v, d.io.vtk_base_pvd, i); continue; }
        if (v[i] == "io.vtk_ext_pvd") { file_read_label_var(v, d.io.vtk_ext_pvd, i); continue; }
        if (v[i] == "io.precision") { file_read_label_var(v, d.io.precision, i); continue; }
        
        // particle information
        if (v[i] == "constant.delta_s") { file_read_label_var(v, d.constant.delta_s, i); continue; }
        
        // fluid distribution
        if (v[i] == "fluid")
        {
            ++i;
            if (v[i] == "block_{")
            {
                // 3d block requires 4 points
                //
                // P3 x
                //    | x P4
                //    |/
                //    x-----x
                //   P1    P2
                
                vec3<double> p1, p2, p3, p4, vel, vel_o;
                double rho = 0.0, rho_o = 0.0, mu = 0.0, gamma = 0.0;
                
                while (++i, v[i] != "}_block")
                {
                    if (v[i] == "p1") { file_read_label_vec(v, p1, i); continue; }
                    if (v[i] == "p2") { file_read_label_vec(v, p2, i); continue; }
                    if (v[i] == "p3") { file_read_label_vec(v, p3, i); continue; }
                    if (v[i] == "p4") { file_read_label_vec(v, p4, i); continue; }
                    if (v[i] == "vel") { file_read_label_vec(v, vel, i); continue; }
                    if (v[i] == "vel_o") { file_read_label_vec(v, vel_o, i); continue; }
                    if (v[i] == "rho") { file_read_label_var(v, rho, i); continue; }
                    if (v[i] == "rho_o") { file_read_label_var(v, rho_o, i); continue; }
                    if (v[i] == "mu") { file_read_label_var(v, mu, i); continue; }
                    if (v[i] == "gamma") { file_read_label_var(v, gamma, i); continue; }
                    file_config_read_error(v[i], counter);
                }
                
                vec3<double> unit_p12 = normalise(p2 - p1);     // unit vector p1 to p2
                vec3<double> unit_p13 = normalise(p3 - p1);     // unit vector p1 to p3
                vec3<double> unit_p14 = normalise(p4 - p1);     // unit vector p1 to p4
                long size_p12 = std::round(mag(p2 - p1) / d.constant.delta_s);       // size of p1 to p2
                long size_p13 = std::round(mag(p3 - p1) / d.constant.delta_s);       // size of p1 to p3
                long size_p14 = std::round(mag(p4 - p1) / d.constant.delta_s);       // size of p1 to p3
                
                for (long ii = 0; ii < size_p12; ++ii)
                {
                    vec3<double> pos1 = p1 + unit_p12 * ii * d.constant.delta_s + unit_p12 * 0.5 * d.constant.delta_s;
                    for (long jj = 0; jj < size_p13; ++jj)
                    {
                        vec3<double> pos2 = pos1 + unit_p13 * jj * d.constant.delta_s + unit_p13 * 0.5 * d.constant.delta_s;
                        for (long kk = 0; kk < size_p14; ++kk)
                        {
                            d.fluid.pos.emplace_back(pos2 + unit_p14 * kk * d.constant.delta_s + unit_p14 * 0.5 * d.constant.delta_s);  // pos at n
                            // calculate position of s_n-1 taking initial velocities into consideration
                            // using s_n-1 = s_n - u_n-1/2 dt where u_n-1/2 = (u_n + u_n-1) / 2
                            // s_n-1 = s_n - dt_n-1 (u_n + u_n-1) / 2
                            d.fluid.pos_o.emplace_back(d.fluid.pos.back() - d.time.dt_max * 0.5 * (vel + vel_o));  // pos at n-1
                            d.fluid.m.emplace_back(rho * std::pow(d.constant.delta_s, d.constant.dim));                   // mass = rho * V
                            d.fluid.vel.emplace_back(vel);                                            // vel at n
                            d.fluid.vel_o.emplace_back(vel_o);                                        // vel at n-1
                            d.fluid.rho.emplace_back(rho);                                            // rho at n
                            d.fluid.rho_o.emplace_back(rho_o);                                        // rho at n-1
                            d.fluid.mu.emplace_back(mu);                                              // mu at n
                            d.fluid.gamma.emplace_back(gamma);                                              // gamma at n
                        }
                    }
                }
            }
            if (v[i] == "plane_{")
            {
                // 3d plane requires 3 points
                //
                //      x P3
                //     /
                //    x-----x
                //   P1    P2
                
                vec3<double> p1, p2, p3, vel, vel_o;
                double rho = 0.0, rho_o = 0.0, mu = 0.0, gamma = 0.0;
                
                while (++i, v[i] != "}_plane")
                {
                    if (v[i] == "p1") { file_read_label_vec(v, p1, i); continue; }
                    if (v[i] == "p2") { file_read_label_vec(v, p2, i); continue; }
                    if (v[i] == "p3") { file_read_label_vec(v, p3, i); continue; }
                    if (v[i] == "vel") { file_read_label_vec(v, vel, i); continue; }
                    if (v[i] == "vel_o") { file_read_label_vec(v, vel_o, i); continue; }
                    if (v[i] == "rho") { file_read_label_var(v, rho, i); continue; }
                    if (v[i] == "rho_o") { file_read_label_var(v, rho_o, i); continue; }
                    if (v[i] == "mu") { file_read_label_var(v, mu, i); continue; }
                    if (v[i] == "gamma") { file_read_label_var(v, gamma, i); continue; }
                    file_config_read_error(v[i], counter);
                }
                
                vec3<double> unit_p12 = normalise(p2 - p1);     // unit vector p1 to p2
                vec3<double> unit_p13 = normalise(p3 - p1);     // unit vector p1 to p3
                long size_p12 = std::round(mag(p2 - p1) / d.constant.delta_s);       // size of p1 to p2
                long size_p13 = std::round(mag(p3 - p1) / d.constant.delta_s);       // size of p1 to p3
                
                for (long ii = 0; ii < size_p12; ++ii)
                {
                    vec3<double> pos1 = p1 + unit_p12 * ii * d.constant.delta_s + unit_p12 * 0.5 * d.constant.delta_s;
                    for (long jj = 0; jj < size_p13; ++jj)
                    {
                        d.fluid.pos.emplace_back(pos1 + unit_p13 * jj * d.constant.delta_s + unit_p13 * 0.5 * d.constant.delta_s);      // pos at n
                        // calculate position of s_n-1 taking initial velocities into consideration
                        // using s_n-1 = s_n - u_n-1/2 dt where u_n-1/2 = (u_n + u_n-1) / 2
                        // s_n-1 = s_n - dt_n-1 (u_n + u_n-1) / 2
                        d.fluid.pos_o.emplace_back(d.fluid.pos.back() - d.time.dt_max * 0.5 * (vel + vel_o));  // pos at n-1
                        d.fluid.m.emplace_back(rho * std::pow(d.constant.delta_s, d.constant.dim));                   // mass = rho * V
                        d.fluid.vel.emplace_back(vel);                                            // vel at n
                        d.fluid.vel_o.emplace_back(vel_o);                                        // vel at n-1
                        d.fluid.rho.emplace_back(rho);                                            // rho at n
                        d.fluid.rho_o.emplace_back(rho_o);                                        // rho at n-1
                        d.fluid.mu.emplace_back(mu);                                              // mu at n
                        d.fluid.gamma.emplace_back(gamma);                                              // mu at n
                    }
                }
            }
            continue;
        }
        
        // boundary distribution
        if (v[i] == "boundary")
        {
            ++i;
            if (v[i] == "plane_{")
            {
                // 3d plane requires 3 points
                //
                //      x P3
                //     /
                //    x-----x
                //    P1    P2
                
                vec3<double> p1, p2, p3, vel;
                while (++i, v[i] != "}_plane")
                {
                    if (v[i] == "p1") { file_read_label_vec(v, p1, i); continue; }
                    if (v[i] == "p2") { file_read_label_vec(v, p2, i); continue; }
                    if (v[i] == "p3") { file_read_label_vec(v, p3, i); continue; }
                    if (v[i] == "vel") { file_read_label_vec(v, vel, i); continue; }
                    file_config_read_error(v[i], counter);
                }
                
                vec3<double> unit_p12 = normalise(p2 - p1);     // unit vector p1 to p2
                vec3<double> unit_p13 = normalise(p3 - p1);     // unit vector p1 to p3
                double delta = 0.5 * d.constant.delta_s;  // 0.5 factor because of the local symmetry boundary condition
                long size_p12 = std::round(mag(p2 - p1) / delta);       // size of p1 to p2
                long size_p13 = std::round(mag(p3 - p1) / delta);       // size of p1 to p3
                for (long ii = 0; ii < size_p12; ++ii)
                {
                    vec3<double> pos = p1 + unit_p12 * ii * delta + unit_p12 * 0.5 * delta;
                    for (long jj = 0; jj < size_p13; ++jj)
                    {
                        d.boundary.pos.emplace_back(pos + unit_p13 * jj * delta + unit_p13 * 0.5 * delta);  // pos_s at n
                        d.boundary.vel.emplace_back(vel);                                                       // vel at n
                    }
                }
            }
            if (v[i] == "line_{")
            {
                // 3d line requires 2 points
                //
                //    x-----x
                //    P1    P2
                
                vec3<double> p1, p2, vel;
                while (++i, v[i] != "}_line")
                {
                    if (v[i] == "p1") { file_read_label_vec(v, p1, i); continue; }
                    if (v[i] == "p2") { file_read_label_vec(v, p2, i); continue; }
                    if (v[i] == "vel") { file_read_label_vec(v, vel, i); continue; }
                    file_config_read_error(v[i], counter);
                }
                
                vec3<double> unit_p12 = normalise(p2 - p1);     // unit vector p1 to p2
                double delta = 0.5 * d.constant.delta_s;  // 0.5 factor because of the local symmetry boundary condition
                long size_p12 = std::round(mag(p2 - p1) / delta);       // size of p1 to p2
                for (long ii = 0; ii < size_p12; ++ii)
                {
                    d.boundary.pos.emplace_back(p1 + unit_p12 * ii * delta + unit_p12 * 0.5 * delta);   // pos_s at n
                    d.boundary.vel.emplace_back(vel);                                                       // vel at n
                }
            }
            continue;
        }
        file_config_read_error(v[i], counter);
    }
}

void file_config_read_error(const std::string &str, long counter)
{
    std::cerr << "Error! Line: " << counter << ", option \"" << str << "\" does not exist." << std::endl;
    exit(EXIT_FAILURE);
}

// generate filename for datafile
const std::string file_core_generate(data &d)
{
    return d.io.core_base + "-" + to_string_with_padding(d.io.itr, '0') + d.io.core_ext;      // generate filename with time info
}

// write datafile
void file_core_write(data &d, const std::string &filename)
{
    filename_rename_ifexist(filename);
    
    std::ofstream ofs;
    
    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (ofs.good())
    {
        ofs.precision(d.io.precision);
        
        file_write_array1dvec(ofs, "fluid.pos", d.fluid.pos);          // array of particle position vectors at n
        file_write_array1dvec(ofs, "fluid.pos_o", d.fluid.pos_o);          // array of particle position vectors at n
        file_write_array1dvec(ofs, "fluid.vel", d.fluid.vel);          // array of particle velocity vectors at n
        file_write_array1dvec(ofs, "fluid.vel_o", d.fluid.vel_o);          // array of particle velocity vectors at n
        file_write_array1dvec(ofs, "boundary.pos", d.boundary.pos);          // array of wall particle position vectors at n
        file_write_array1dvec(ofs, "boundary.vel", d.boundary.vel);          // array of wall particle velocity vectors at n
        
        file_write_array1d(ofs, "fluid.rho", d.fluid.rho);          // array of particle density at n
        file_write_array1d(ofs, "fluid.rho_o", d.fluid.rho_o);          // array of particle density at n
        file_write_array1d(ofs, "fluid.m", d.fluid.m);              // array of particle mass
        file_write_array1d(ofs, "fluid.mu", d.fluid.mu);              // array of particle viscosity
        file_write_array1d(ofs, "fluid.gamma", d.fluid.gamma);              // array of particle surface tension
        
        file_write_vec(ofs, "constant.g", d.constant.g);                          // gravity
        
        file_write_var(ofs, "constant.h_eta", d.constant.h_eta);                   // expansion ratio for smoothing length to capture sufficient neighbouring particles
        file_write_var(ofs, "constant.h_kappa", d.constant.h_kappa);               // support domain to smoothing length ratio
        file_write_var(ofs, "constant.h", d.constant.h);                       // constant smoothing length for each particle
        
        file_write_var(ofs, "time.t_begin", d.time.t_begin);           // default start time
        file_write_var(ofs, "time.t_end", d.time.t_end);               // default end time
        file_write_var(ofs, "time.t", d.time.t);                       // current time
        file_write_var(ofs, "time.dt", d.time.dt);                     // time step
        file_write_var(ofs, "time.dt_o", d.time.dt_o);                     // time step
        file_write_var(ofs, "time.dt_c_cfl", d.time.dt_c_cfl);         // time step cfl coefficient
        file_write_var(ofs, "time.dt_c_viscosity", d.time.dt_c_viscosity);         // time step viscosity coefficient
        file_write_var(ofs, "time.dt_c_surfacetension", d.time.dt_c_surfacetension);         // time step surface tension coefficient
        file_write_var(ofs, "time.dt_max", d.time.dt_max);             // max time step
        file_write_var(ofs, "time.dt_output", d.time.dt_output);         // time for next output
        
        file_write_var(ofs, "constant.rho_0", d.constant.rho_0);                      // reference density
        file_write_var(ofs, "constant.c_0", d.constant.c_0);                          // sound speed factor
        
        file_write_var(ofs, "constant.delta_s", d.constant.delta_s);                  // default particle spacing
        
        file_write_var(ofs, "io.core_base", d.io.core_base);          // datafile basename
        file_write_var(ofs, "io.core_ext", d.io.core_ext);            // datafile extension
        file_write_var(ofs, "io.vtk_base_vtu", d.io.vtk_base_vtu);                    // vtu basename
        file_write_var(ofs, "io.vtk_ext_vtu", d.io.vtk_ext_vtu);                      // vtu extension
        file_write_var(ofs, "io.vtk_base_pvd", d.io.vtk_base_pvd);                    // vtu basename
        file_write_var(ofs, "io.vtk_ext_pvd", d.io.vtk_ext_pvd);                      // vtu extension
        file_write_var(ofs, "kernel.type", d.kernel.type);                        // kernel type
        
        file_write_var(ofs, "constant.dim", d.constant.dim);                // dimensional number
        file_write_var(ofs, "io.itr", d.io.itr);                // output iterator number
        
        file_write_var(ofs, "io.precision", d.io.precision);            // input/output precision
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    ofs.close();
}

// read datafile
void file_core_read(data &d, const std::string &filename)
{
    std::ifstream ifs;
    
    ifs.open(filename);
    if (ifs.is_open())                                  // Check if file is open
    {
        std::string line;
        while (std::getline(ifs, line))
        {
            file_core_read_parse(d, line);
        }
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
}

// read parse lines in datafile
void file_core_read_parse(data &d, const std::string &line)
{
    std::vector<std::string> v = tokenise<std::string>(line, " \t");    // tokenise using space and tab
    for (long i = 0, i_size=v.size(); i != i_size; ++i)
    {                                                                   // check token is
        if ((v[i] == "#") ||                                            // #
            (v[i] == ";") ||                                            // ;
            (v.empty()) ||                                              // empty line
            (v[i][0] == '#') ||                                         // begin with #
            (v[i][0] == ';'))                                           // begin with ;
        {
            v.erase(v.begin() + i, v.begin() + i_size);                 // erase token and everthing after it
            break;                                                      // break from for loop
        }
        std::transform(v[i].begin(), v[i].end(), v[i].begin(), ::tolower); // converts to lower case
    }
    file_core_read_label(d, v);
}

// read options
void file_core_read_label(data &d, const std::vector<std::string> &v)
{
    for (long i = 0, i_size = v.size(); i < i_size; ++i)
    {
        if (v[i] == "fluid.pos") { file_read_label_array1dvec(v, d.fluid.pos, i); continue; }
        if (v[i] == "fluid.pos_o") { file_read_label_array1dvec(v, d.fluid.pos_o, i); continue; }
        if (v[i] == "fluid.vel") { file_read_label_array1dvec(v, d.fluid.vel, i); continue; }
        if (v[i] == "fluid.vel_o") { file_read_label_array1dvec(v, d.fluid.vel_o, i); continue; }
        if (v[i] == "boundary.pos") { file_read_label_array1dvec(v, d.boundary.pos, i); continue; }
        if (v[i] == "boundary.vel") { file_read_label_array1dvec(v, d.boundary.vel, i); continue; }
        
        if (v[i] == "fluid.rho") { file_read_label_array1d(v, d.fluid.rho, i); continue; }
        if (v[i] == "fluid.rho_o") { file_read_label_array1d(v, d.fluid.rho_o, i); continue; }
        if (v[i] == "fluid.m") { file_read_label_array1d(v, d.fluid.m, i); continue; }
        if (v[i] == "fluid.mu") { file_read_label_array1d(v, d.fluid.mu, i); continue; }
        if (v[i] == "fluid.gamma") { file_read_label_array1d(v, d.fluid.gamma, i); continue; }
        
        if (v[i] == "constant.g") { file_read_label_vec(v, d.constant.g, i); continue; }
        
        if (v[i] == "constant.h_eta") { file_read_label_var(v, d.constant.h_eta, i); continue; }
        if (v[i] == "constant.h_kappa") { file_read_label_var(v, d.constant.h_kappa, i); continue; }
        if (v[i] == "constant.h") { file_read_label_var(v, d.constant.h, i); continue; }
        
        if (v[i] == "time.t_begin") { file_read_label_var(v, d.time.t_begin, i); continue; }
        if (v[i] == "time.t_end") { file_read_label_var(v, d.time.t_end, i); continue; }
        if (v[i] == "time.t") { file_read_label_var(v, d.time.t, i); continue; }
        if (v[i] == "time.dt") { file_read_label_var(v, d.time.dt, i); continue; }
        if (v[i] == "time.dt_o") { file_read_label_var(v, d.time.dt_o, i); continue; }
        if (v[i] == "time.dt_c_cfl") { file_read_label_var(v, d.time.dt_c_cfl, i); continue; }
        if (v[i] == "time.dt_c_viscosity") { file_read_label_var(v, d.time.dt_c_viscosity, i); continue; }
        if (v[i] == "time.dt_c_surfacetension") { file_read_label_var(v, d.time.dt_c_surfacetension, i); continue; }
        if (v[i] == "time.dt_max") { file_read_label_var(v, d.time.dt_max, i); continue; }
        if (v[i] == "time.dt_output") { file_read_label_var(v, d.time.dt_output, i); continue; }
        
        if (v[i] == "constant.rho_0") { file_read_label_var(v, d.constant.rho_0, i); continue; }
        if (v[i] == "constant.c_0") { file_read_label_var(v, d.constant.c_0, i); continue; }
        
        if (v[i] == "constant.delta_s") { file_read_label_var(v, d.constant.delta_s, i); continue; }
        
        if (v[i] == "io.core_base") { file_read_label_var(v, d.io.core_base, i); continue; }
        if (v[i] == "io.core_ext") { file_read_label_var(v, d.io.core_ext, i); continue; }
        if (v[i] == "io.vtk_base_vtu") { file_read_label_var(v, d.io.vtk_base_vtu, i); continue; }
        if (v[i] == "io.vtk_ext_vtu") { file_read_label_var(v, d.io.vtk_ext_vtu, i); continue; }
        if (v[i] == "io.vtk_base_pvd") { file_read_label_var(v, d.io.vtk_base_pvd, i); continue; }
        if (v[i] == "io.vtk_ext_pvd") { file_read_label_var(v, d.io.vtk_ext_pvd, i); continue; }
        if (v[i] == "kernel.type") { file_read_label_var(v, d.kernel.type, i); continue; }
        
        if (v[i] == "constant.dim") { file_read_label_var(v, d.constant.dim, i); continue; }
        if (v[i] == "io.itr") { file_read_label_var(v, d.io.itr, i); continue; }
        
        if (v[i] == "io.precision") { file_read_label_var(v, d.io.precision, i); continue; }
        
        file_core_read_error(v[i]);
    }
}

void file_core_read_error(const std::string &str)
{
    std::cerr << "Error! reading file. Option " << str << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
}

const std::string file_vtk_vtu_generate(data &d)
{
    return d.io.vtk_base_vtu + "-" + to_string_with_padding(d.io.itr, '0') + d.io.vtk_ext_vtu;      // generate vtu filename with time info
}

const std::string file_vtk_pvd_generate(data &d)
{
    return d.io.vtk_base_pvd + d.io.vtk_ext_pvd;      // generate vtu filename with time info
}

void file_vtk_vtu_write(data &d, const std::string &filename)
{
    filename_rename_ifexist(filename);                          // rename filename if exist
    
    std::ofstream ofs;
    
    ofs.open(filename, std::ofstream::trunc);
    
    if (ofs.good())
    {
        ofs.precision(d.io.precision);                                      // set output precision
        long array_size = d.fluid.pos.size() + d.boundary.pos.size();           // array_size = fluid + wall particles
        auto dim = d.fluid.pos[0].data.size();                              // get dimensions
        
        // file header
        ofs << "<?xml version=\"1.0\"?>" << std::endl;
        
        // open vtkfile
        ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        
        // open unstructuredgrid
        ofs << "  <UnstructuredGrid>" << std::endl;
        
        // open piece
        ofs << "    <Piece NumberOfPoints=\"" << array_size << "\" NumberOfCells=\"" << array_size << "\">" << std::endl;
        
        // open pointdata and specify default active pointdata
        ofs << "      <PointData Scalars=\"Index\" Vectors=\"Velocity\">" << std::endl;
        
        // pressure
        ofs << "        <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (auto i: d.fluid.p)
        {
            ofs << i << " ";
        }
        for (long i = 0, i_size = d.boundary.pos.size(); i < i_size; ++i)
        {
            ofs << "0" << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // density
        ofs << "        <DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (auto i: d.fluid.rho)
        {
            ofs << i << " ";
        }
        for (long i = 0, i_size = d.boundary.pos.size(); i < i_size; ++i)
        {
            ofs << "0" << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // velocity
        ofs << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (auto i: d.fluid.vel)
        {
            ofs << i << " ";
        }
        for (long i = 0, i_size = d.boundary.pos.size(); i < i_size; ++i)
        {
            ofs << "0 0 0" << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // acceleration
        ofs << "        <DataArray type=\"Float32\" Name=\"Acceleration\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (auto i: d.fluid.acc)
        {
            ofs << i << " ";
        }
        for (long i = 0, i_size = d.boundary.pos.size(); i < i_size; ++i)
        {
            ofs << "0 0 0" << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
                
        // close pointdata
        ofs << "      </PointData>" << std::endl;
        
        // open points
        ofs << "      <Points>" << std::endl;
        
        // position
        ofs << "        <DataArray type=\"Float32\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (auto i: d.fluid.pos)
        {
            ofs << i << " ";
        }
        for (auto i: d.boundary.pos)
        {
            ofs << i << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // close points
        ofs << "      </Points>" << std::endl;
        
        // open cells
        ofs << "      <Cells>" << std::endl;
        
        // connectivity
        ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (long i = 0; i != array_size; ++i)
        {
            ofs << i << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // offsets
        ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (long i = 0; i != array_size; ++i)
        {
            ofs << i + 1 << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // types
        ofs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
        ofs << "          ";
        for (long i = 0; i != array_size; ++i)
        {
            ofs << "1" << " ";
        }
        ofs << std::endl;
        ofs << "        </DataArray>" << std::endl;
        
        // close cells
        ofs << "      </Cells>" << std::endl;
        
        // close piece
        ofs << "    </Piece>" << std::endl;
        
        // close unstructuredgrid
        ofs << "  </UnstructuredGrid>" << std::endl;
        
        // close vtkfile
        ofs << "</VTKFile>" << std::endl;
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs.close();
}

void file_vtk_pvd_write(data &d, const std::string &pvdfilename, const std::string &vtufilename)
{
    if (d.io.itr == 0)                                  // new simulation
    {
        filename_rename_ifexist(pvdfilename);               // rename filename if exist
        file_vtk_pvd_write_begin(pvdfilename);                   // create new file
    }
    
    file_vtk_pvd_write_app(d, pvdfilename, vtufilename);         // append to existing file
    
    if (near_equal(d.time.t, d.time.t_end))         // end of simulation
    {
        file_vtk_pvd_write_end(pvdfilename);                     // create new file
    }
}

void file_vtk_pvd_write_begin(const std::string &filename)
{
    std::ofstream ofs;
    
    ofs.open(filename, std::ofstream::trunc);
    
    if (ofs.good())
    {
        // file header
        ofs << "<?xml version=\"1.0\"?>" << std::endl;
        
        // open vtkfile
        ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        
        // open collection
        ofs << "  <Collection>" << std::endl;
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs.close();
}

void file_vtk_pvd_write_app(data &d, const std::string &pvdfilename, const std::string &vtufilename)
{
    std::ofstream ofs;
    
    ofs.open(pvdfilename, std::ofstream::app);
    
    if (ofs.good())
    {
        // insert dataset info
        ofs << "    <DataSet timestep=\"" << d.time.t << "\" group=\"\" part=\"0\" file=\"" << vtufilename << "\"/>" << std::endl;
    }
    else
    {
        std::cerr << pvdfilename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs.close();
}

void file_vtk_pvd_write_end(const std::string &filename)
{
    std::ofstream ofs;
    
    ofs.open(filename, std::ofstream::app);
    
    if (ofs.good())
    {
        // close collection
        ofs << "  </Collection>" << std::endl;
        
        // close vtkfile
        ofs << "</VTKFile>" << std::endl;
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs.close();
}

void file_vtk_pvd_reformat(data &d)
{
    std::string filename = file_vtk_pvd_generate(d);
    std::ifstream ifs;
    std::ofstream ofs;
    std::string temp_filename = filename + ".tmp";
    
    ifs.open(filename);
    ofs.open(temp_filename, std::ofstream::trunc);
    
    if (ifs.is_open())                                  // Check if file is open
    {
        std::string line;
        std::string target = "timestep=\"" + to_string(d.time.t) + "\"" ;           // target keyword to stop reading
        
        while (getline(ifs, line))                      // read line by line
        {
            if (line.find(target) != std::string::npos) // if target found
            {
                break;                                  // exit while loop
            }
            
            if (ofs.good())
            {
                ofs << line << std::endl;               // copy line to temp file
            }
            else
            {
                std::cerr << temp_filename << " error: " << std::strerror(errno) << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs.close();
    ifs.close();
    
    if (std::rename(temp_filename.c_str(), filename.c_str()))           // rename temp file to original file
    {
        std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
}
