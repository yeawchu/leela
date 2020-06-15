//
//  output.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "output.hpp"

void output(data &d)
{
    output_file(d);
    output_screen(d);
}

void output_file(data &d)
{
    if (near_equal(d.time.t, static_cast<double>(d.io.itr) * d.time.dt_output))
    {
        std::string datafile_name = file_core_generate(d);   // generate input filename
        file_core_write(d, datafile_name);                   // write to file
        
        std::string vtu_name = file_vtk_vtu_generate(d);             // generate input filename;
        file_vtk_vtu_write(d, vtu_name);                             // write to file
        
        std::string pvd_name = file_vtk_pvd_generate(d);             // generate input filename;
        file_vtk_pvd_write(d, pvd_name, vtu_name);                             // write to file
    }
}

void output_screen(data &d)
{
    output_screen_progress(d);
}

void output_screen_progress(data &d)
{
    output_screen_progress_minimal(d);
}

void output_screen_progress_bar(data &d)
{
    double fraction = (d.time.t - d.time.t_begin) / d.time.t_end;
    int bar_width = 70;
    
    std::cout << "[";
    int pos = bar_width * fraction;
    for (int i = 0; i != bar_width; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    int precision = 6;
    std::cout << "] " << std::fixed << std::setprecision(precision) << (fraction * 100.0) << " %\r";
    std::cout.flush();
}

void output_screen_progress_info(data &d)
{
    double fraction = (d.time.t - d.time.t_begin) / d.time.t_end;
    int precision = 9;
    std::cout << std::fixed << std::setprecision(precision) << "[ t " << d.time.t << " ] [ dt " << d.time.dt << " ] [ % " << (fraction * 100.0) << " ]\r";
    std::cout.flush();
}

void output_screen_progress_minimal(data &d)
{
    if (near_equal(d.time.t / d.time.dt_output, static_cast<double>(d.io.itr)))
    {
        double fraction = (d.time.t - d.time.t_begin) / d.time.t_end;
        int precision = 9;
        std::cout << "[ io " << d.io.itr << " ] " << std::fixed << std::setprecision(precision) << "[ t " << d.time.t << " ] [ % " << (fraction * 100.0) << " ]\r";
        std::cout.flush();
    }
}

void output_screen_progress_debug(data &d)
{
    double fraction = (d.time.t - d.time.t_begin) / d.time.t_end;
    int precision = 9;
    std::cout << std::fixed << std::setprecision(precision) << "[ t " << d.time.t << " ] [ dt " << d.time.dt << " ] [ % " << (fraction * 100.0) << " ]" << std::endl;
}

