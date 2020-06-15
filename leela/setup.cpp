//
//  setup.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "setup.hpp"

void setup(data &d)
{
    setup_neighbour(d);
    setup_kernel(d);
}

void setup_neighbour(data &d)
{
    setup_neighbour_search(d);
}

void setup_neighbour_search(data &d)
{
    //setup_neighbour_search_cell(d);
    setup_neighbour_search_cellpair(d);
}

void setup_neighbour_search_cell(data &d)
{
    // Create fixed smoothing length cell search,
    //      |---|---|---|
    //      |   |   |   |
    //      |---|---|---|
    //      |   | x |   |
    //      |---|---|---|
    //      |   |   |   |
    //      |---|---|---|

    std::vector<vec3<double>> particle_pos;     // holds all particle positions
    particle_pos.reserve(d.fluid.pos.size() + d.boundary.pos.size());
    particle_pos.insert(particle_pos.end(), d.fluid.pos.begin(), d.fluid.pos.end());        // insert fluid particles
    particle_pos.insert(particle_pos.end(), d.boundary.pos.begin(), d.boundary.pos.end());    // insert solid particles
    
    vec3<double> domain_lb = particle_pos[0];    // set initial domain lower and upper bounds to reference particle position
    vec3<double> domain_ub = particle_pos[0];
    
    for (auto &pos: particle_pos)      // determine domain lower and upper bound limits
    {
        domain_lb = min(domain_lb, pos);
        domain_ub = max(domain_ub, pos);
    }

    vec3<double> particle_volume(d.constant.delta_s, d.constant.delta_s, d.constant.delta_s);    // set particle volume
    
    vec3<double> cell_lb = domain_lb - 0.5 * particle_volume;   // set cell lower and upper bound limits
    vec3<double> cell_ub = domain_ub + 0.5 * particle_volume;
    vec3<double> cell_bound_length = cell_ub - cell_lb;     // cell length between bound limits
    double cell_length = d.constant.h_kappa * d.constant.h;       // set cell length = h_kappa * h
    vec3<long> cell_size = cast<long>(ceil(cell_bound_length / cell_length));   // work out number of cells in each direction
    
    std::vector<std::vector<long>> cell_particle_index;
    cell_particle_index.resize(cell_size.x * cell_size.y * cell_size.z);       // resize outer column of 2D array
    long particle_index_size = d.constant.dim * static_cast<long>(cell_length / d.constant.delta_s);    // estimate number of particles in a cell

    for (long i = 0, i_size = cell_particle_index.size(); i < i_size; ++i)
    {
        cell_particle_index[i].reserve(particle_index_size);       // reserve memory
    }

    for (long i = 0, i_size = particle_pos.size(); i != i_size; ++i)       // for each particle position:
    {
        vec3<long> particle_cell_pos = cast<long>(floor((particle_pos[i] - cell_lb) / cell_length));         // determine cell position
        long cell_index = cell_pos_to_idx<long>(particle_cell_pos, cell_size);               // find cell index
        cell_particle_index[cell_index].emplace_back(i);            // store particle index in cell
    }
    
    d.nlist.ff.i.clear();                          // clear fluid to fluid particle neighbour list
    d.nlist.ff.j.clear();                          // clear fluid to fluid particle neighbour list
    d.nlist.fb.i.clear();                          // clear fluid to wall particle neighbour list
    d.nlist.fb.j.clear();                          // clear fluid to wall particle neighbour list
    long fluid_size = d.fluid.pos.size();           // get number of fluid particles
    
    long nl_size = static_cast<long>(particle_index_size * fluid_size);        // estimate of nl size
    d.nlist.ff.i.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.ff.j.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.fb.i.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.fb.j.reserve(nl_size);     // need to reserve d.nl.ff

    for (long ci = 0, ci_size = cell_particle_index.size(); ci != ci_size; ++ci)           // for each cell:
    {
        std::vector<long> c_nlist = cell_nlist<long>(ci, cell_size);   // determine adjacent cell neighbour list

        for (auto &cni: c_nlist)                    // for each adjacent cell neighbour:
        {
            for (long i = 0, i_size = cell_particle_index[ci].size(); i != i_size; ++i)    // for each particle in cell:
            {
                for (long j = 0, j_size = cell_particle_index[cni].size(); j != j_size; ++j)      // for each neighbour particle index in cell neighbour list:
                {
                    long pi = cell_particle_index[ci][i];                // get particle index in cell
                    long pni = cell_particle_index[cni][j];      // get neighbour particle index in the same or neighbouring cells
                    
                    if (dist2(particle_pos[pi], particle_pos[pni]) <= cell_length * cell_length)  // is particle neighbour within range?
                    {
                        if (pi < fluid_size)               // if particle is fluid
                        {
                            if (pni < fluid_size)           // if neighbour particle is fluid
                            {
                                d.nlist.ff.i.emplace_back(pi);
                                d.nlist.ff.j.emplace_back(pni);
                            }
                            else                            // if neighbour particle is not fluid
                            {
                                d.nlist.fb.i.emplace_back(pi);
                                d.nlist.fb.j.emplace_back(pni - fluid_size);
                            }
                        }
                    }
                }
            }
        }
    }
}

void setup_neighbour_search_cellpair(data &d)
{
    // Performs fixed smoothing length cell search using pair approach,
    // ie only search in the forward direction (see Fig)
    //      |---|---|---|
    //      |   |   |   |
    //      |---|---|---|
    //          | x |   |
    //          |---|---|

    std::vector<vec3<double>> particle_pos;     // holds all particle positions
    particle_pos.reserve(d.fluid.pos.size() + d.boundary.pos.size());
    particle_pos.insert(particle_pos.end(), d.fluid.pos.begin(), d.fluid.pos.end());        // insert fluid particles
    particle_pos.insert(particle_pos.end(), d.boundary.pos.begin(), d.boundary.pos.end());    // insert solid particles

    vec3<double> domain_lb = particle_pos[0];    // set initial domain lower and upper bounds to reference particle position
    vec3<double> domain_ub = particle_pos[0];
    
    for (auto &pos: particle_pos)      // determine domain lower and upper bound limits
    {
        domain_lb = min(domain_lb, pos);
        domain_ub = max(domain_ub, pos);
    }
    vec3<double> particle_volume(d.constant.delta_s, d.constant.delta_s, d.constant.delta_s);    // set particle volume
    
    vec3<double> cell_lb = domain_lb - 0.5 * particle_volume;   // set cell lower and upper bound limits
    vec3<double> cell_ub = domain_ub + 0.5 * particle_volume;
    vec3<double> cell_bound_length = cell_ub - cell_lb;     // cell length between bound limits
    double cell_length = d.constant.h_kappa * d.constant.h;       // set cell length = h_kappa * h
    vec3<long> cell_size = cast<long>(ceil(cell_bound_length / cell_length));   // work out number of cells in each direction

    std::vector<std::vector<long>> cell_particle_index;
    cell_particle_index.resize(cell_size.x * cell_size.y * cell_size.z);       // resize outer column of 2D array
    long particle_index_size = d.constant.dim * static_cast<long>(cell_length / d.constant.delta_s);    // estimate number of particles in a cell

    for (long i = 0, i_size = cell_particle_index.size(); i < i_size; ++i)
    {
        cell_particle_index[i].reserve(particle_index_size);       // reserve memory
    }

    for (long i = 0, i_size = particle_pos.size(); i != i_size; ++i)       // for each particle position:
    {
        vec3<long> particle_cell_pos = cast<long>(floor((particle_pos[i] - cell_lb) / cell_length));         // determine cell position
        long cell_index = cell_pos_to_idx<long>(particle_cell_pos, cell_size);               // find cell index
        cell_particle_index[cell_index].emplace_back(i);            // store particle index in cell
    }

    d.nlist.ff.i.clear();                          // clear fluid to fluid particle neighbour list
    d.nlist.ff.j.clear();                          // clear fluid to fluid particle neighbour list
    d.nlist.fb.i.clear();                          // clear fluid to wall particle neighbour list
    d.nlist.fb.j.clear();                          // clear fluid to wall particle neighbour list
    long fluid_size = d.fluid.pos.size();           // get number of fluid particles
    
    long nl_size = static_cast<long>(particle_index_size * fluid_size);        // estimate of nl size
    d.nlist.ff.i.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.ff.j.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.fb.i.reserve(nl_size);     // need to reserve d.nl.ff
    d.nlist.fb.j.reserve(nl_size);     // need to reserve d.nl.ff

    for (long ci = 0, ci_size = cell_particle_index.size(); ci != ci_size; ++ci)           // for each cell:
    {
        std::vector<long> c_nlist = cellpair_nlist<long>(ci, cell_size);   // determine adjacent cell neighbour list
        for (auto &cni: c_nlist)                    // for each adjacent cell neighbour:
        {
            for (long i = 0, i_size = cell_particle_index[ci].size(); i != i_size; ++i)    // for each particle in cell:
            {
                long jstart = 0;                        // index of neighbour particle
                if (ci == cni) jstart = i;              // if particle is in the same cell, set particle neighbour index to current particle index to avoid unnecessary repeating/redoing existing paired searches when performing the forward pair-wise search

                for (long j = jstart, j_size = cell_particle_index[cni].size(); j != j_size; ++j)      // for each neighbour particle index in cell neighbour list:
                {
                    long pi = cell_particle_index[ci][i];                // get particle index in cell
                    long pni = cell_particle_index[cni][j];      // get neighbour particle index in the same or neighbouring cells
                    
                    if (dist2(particle_pos[pi], particle_pos[pni]) <= cell_length * cell_length)  // is particle neighbour within range?
                    {
                        if (pi < fluid_size)               // if particle is fluid
                        {
                            if (pni < fluid_size)           // if neighbour particle is fluid
                            {
                                d.nlist.ff.i.emplace_back(pi);
                                d.nlist.ff.j.emplace_back(pni);
                                if (pi == pni) continue;                            // prevent from copying itself twice
                                d.nlist.ff.i.emplace_back(pni);
                                d.nlist.ff.j.emplace_back(pi);
                            }
                            else                            // if neighbour particle is not fluid
                            {
                                d.nlist.fb.i.emplace_back(pi);
                                d.nlist.fb.j.emplace_back(pni - fluid_size);
                            }
                        }
                        else                                // if particle is wall
                        {
                            if (pni < fluid_size)
                            {
                                d.nlist.fb.i.emplace_back(pni);
                                d.nlist.fb.j.emplace_back(pi - fluid_size);
                            }
                        }
                    }
                }
            }
        }
    }
}

void setup_kernel(data &d)
{
    d.nlist.ff.w.clear();
    d.nlist.ff.w_dash.clear();
    d.nlist.ff.w_e.clear();
    d.nlist.fb.w.clear();
    d.nlist.fb.w_dash.clear();
    d.nlist.fb.w_e.clear();

    long ff_size = d.nlist.ff.i.size();
    long fb_size = d.nlist.fb.i.size();
    d.nlist.ff.w.resize(ff_size);
    d.nlist.ff.w_dash.resize(ff_size);
    d.nlist.ff.w_e.resize(ff_size);
    d.nlist.fb.w.resize(fb_size);
    d.nlist.fb.w_dash.resize(fb_size);
    d.nlist.fb.w_e.resize(fb_size);

    for (long a = 0, a_size = d.nlist.ff.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.ff.i[a];
        long j = d.nlist.ff.j[a];
        vec3<double> r = d.fluid.pos[i] - d.fluid.pos[j];                               // r_ab
        d.nlist.ff.w[a] = d.kernel.func_w(r, d.kernel.alpha, d.constant.h, d.constant.dim);
        d.nlist.ff.w_dash[a] = d.kernel.func_w_dash(r, d.kernel.alpha, d.constant.h, d.constant.dim);
        d.nlist.ff.w_e[a] = d.kernel.func_w_e(r);
    }

    for (long a = 0, a_size = d.nlist.fb.i.size(); a != a_size; ++a)
    {
        long i = d.nlist.fb.i[a];
        long j = d.nlist.fb.j[a];
        vec3<double> pos_k = 2.0 * d.boundary.pos[j] - d.fluid.pos[i];
        vec3<double> r = d.fluid.pos[i] - pos_k;                  // r_ab using reflective BC, x_k = 2 x_b - x_i
        d.nlist.fb.w[a] = d.kernel.func_w(r, d.kernel.alpha, d.constant.h, d.constant.dim);
        d.nlist.fb.w_dash[a] = d.kernel.func_w_dash(r, d.kernel.alpha, d.constant.h, d.constant.dim);
        d.nlist.fb.w_e[a] = d.kernel.func_w_e(r);
    }
}
