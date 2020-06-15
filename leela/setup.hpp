//
//  setup.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef setup_hpp
#define setup_hpp

#include <iostream>
#include <vector>
#include "vec.h"
#include "data.h"
#include "utils.hpp"

void setup(data &d);
void setup_neighbour(data &d);
void setup_neighbour_search(data &d);
void setup_neighbour_search_cell(data &d);
void setup_neighbour_search_cellpair(data &d);
void setup_kernel(data &d);

//
// Helper template functions for 3D cellpair search
//

// return 1D cell index from vec3 vector
template <typename T>
T cell_pos_to_idx(const vec3<T> &v, const vec3<T> &v_size)
{
    return v.x + (v_size.x * v.y) + (v_size.x * v_size.y * v.z);
}

// return vec3 vector from 1D cell index
template <typename T>
vec3<T> cell_idx_to_pos(const T idx, const vec3<T> &v_size)
{
    T z = static_cast<T>(std::floor(idx / (v_size.x * v_size.y)));
    T y = static_cast<T>(std::floor((idx - z * v_size.x * v_size.y) / v_size.x));
    T x = idx - y * v_size.x - z * v_size.x * v_size.y;
    return vec3<T>(x, y, z);
}

// return downstream neighbouring cells in all directions
template <typename T>
std::vector<T> cell_nlist(const T idx, const vec3<T> &v_size)
{
    vec3<T> v = cell_idx_to_pos(idx, v_size);
    
    T x_min = v.x - 1;
    T y_min = v.y - 1;
    T z_min = v.z - 1;
    T x_max = v.x + 1;
    T y_max = v.y + 1;
    T z_max = v.z + 1;
    
    if (x_min < 0) x_min = 0;
    if (x_max > (v_size.x - 1)) x_max = v_size.x - 1;
    if (y_min < 0) y_min = 0;
    if (y_max > (v_size.y - 1)) y_max = v_size.y - 1;
    if (z_min < 0) z_min = 0;
    if (z_max > (v_size.z - 1)) z_max = v_size.z - 1;
    
    std::vector<T> temp;
    T size = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    temp.reserve(size);
    
    for (T i = x_min; i <= x_max; ++i)
    {
        for (T j = y_min; j <= y_max; ++j)
        {
            for (T k = z_min; k <= z_max; ++k)
            {
                temp.emplace_back(cell_pos_to_idx(vec3<T> (i, j, k), v_size));
            }
        }
    }
    return temp;
}

// return downstream neighbouring cells in all directions
template <typename T>
std::vector<T> cellpair_nlist(const T idx, const vec3<T> &v_size)
{
    vec3<T> v = cell_idx_to_pos(idx, v_size);
    
    T x_min = v.x - 1;
    T y_min = v.y - 1;
    T x_max = v.x + 1;
    T y_max = v.y + 1;
    T z_max = v.z + 1;
    
    if (x_min < 0) x_min = 0;
    if (x_max > (v_size.x - 1)) x_max = v_size.x - 1;
    if (y_min < 0) y_min = 0;
    if (y_max > (v_size.y - 1)) y_max = v_size.y - 1;
    if (z_max > (v_size.z - 1)) z_max = v_size.z - 1;
    
    std::vector<T> temp;
    T size = (x_max - x_min) * (y_max - y_min) * (z_max - v.z);
    temp.reserve(size);
    
    for (T i = x_min; i <= x_max; ++i)
    {
        for (T j = y_min; j <= y_max; ++j)
        {
            for (T k = v.z; k <= z_max; ++k)
            {
                if ((i < v.x) && (j == v.y) && (k == v.z)) continue;
                if ((j < v.y) && (k == v.z)) continue;
                temp.emplace_back(cell_pos_to_idx(vec3<T> (i, j, k), v_size));
            }
        }
    }
    return temp;
}

#endif /* setup_hpp */
