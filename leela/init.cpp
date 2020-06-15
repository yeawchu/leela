//
//  init.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "init.hpp"

// initialise solver
void init(data &d, int argc, const char * argv[])
{
    if (argc == 1)                                  // if no core datafile provided, use config file
    {
        std::string me;                             // define my name
        me = argv[0];                               // get my name
        std::string filename = me + ".config";      // name of config file
        file_config_read(d, filename);              // read init file
        init_config(d);                             // initialise unset file_config parameters
    }
    else if (argc == 2)                             // if core datafile provided
    {
        std::string filename;                       // define filename
        filename = argv[1];                         // get core datafile name
        file_core_read(d, filename);                // read core datafile
        file_vtk_pvd_reformat(d);                   // reformat pvd file for appending
    }
    else                                            // if none of the above conditions are met
    {
        std::string me = argv[0];                   // get my name
        std::cerr << "Usage: \"" << me << "\" or \"" << me << "\" <datafile>" << std::endl;     // return usage message
        exit(EXIT_FAILURE);                         // exit with error code
    }
    init_array(d);                                  // initialise undefined arrays
    init_kernel(d);                                 // initialise kernel parameters and functions
}

// initialise unset config options
void init_config(data &d)
{
    d.constant.h = d.constant.h_eta * d.constant.delta_s;  // smoothing length = h_eta * delta
    
    d.time.t = d.time.t_begin;                                // current time
    d.time.dt = d.time.dt_max;                                // default time step
    d.time.dt_o = d.time.dt_max;                              // default old time step
    
    d.io.itr = 0;                               // start output iterator counter
}

// initialise undefined arrays
void init_array(data &d)
{
    long num = d.fluid.pos.size();                        // get total number of fluid particles
    d.fluid.pos_n.resize(num);                            // define pos_n array size
    d.fluid.vel_n.resize(num);                            // define vel_n array size
    d.fluid.rho_n.resize(num);                            // define rho_n array size
    d.fluid.acc.resize(num);                              // define acc array size
    d.fluid.p.resize(num);                                // define p array size
}

// initialise kernel
void init_kernel(data &d)
{
    if (d.kernel.type == "wendland")
    {
        if (d.constant.dim == 1) d.kernel.alpha = kernel_wendland_alpha_1d();
        if (d.constant.dim == 2) d.kernel.alpha = kernel_wendland_alpha_2d();
        if (d.constant.dim == 3) d.kernel.alpha = kernel_wendland_alpha_3d();
        d.kernel.func_w = kernel_wendland;
        d.kernel.func_w_dash = kernel_wendland_dash;
        d.kernel.func_w_e = kernel_e;
    }
    else
    {
        std::cerr << "W_type: " << d.kernel.type << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}