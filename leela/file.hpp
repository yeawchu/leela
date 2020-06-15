//
//  file.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef file_hpp
#define file_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstring>
#include "data.h"
#include "vec.h"
#include "utils.hpp"

void file_config_read(data &d, const std::string &filename);
void file_config_read_parse(data &d, const std::string &line, long counter);
void file_config_read_label(data &d, const std::vector<std::string> &v, long counter);
void file_config_read_error(const std::string &str, long counter);
const std::string file_core_generate(data &d);
void file_core_write(data &d, const std::string &filename);
void file_core_read(data &d, const std::string &filename);
void file_core_read_parse(data &d, const std::string &line);
void file_core_read_label(data &d, const std::vector<std::string> &v);
void file_core_read_error(const std::string &str);
const std::string file_vtk_vtu_generate(data &d);
const std::string file_vtk_pvd_generate(data &d);
void file_vtk_vtu_write(data &d, const std::string &filename);
void file_vtk_pvd_write(data &d, const std::string &pvdfilename, const std::string &vtufilename);
void file_vtk_pvd_write_begin(const std::string &filename);
void file_vtk_pvd_write_app(data &d, const std::string &pvdfilename, const std::string &vtufilename);
void file_vtk_pvd_write_end(const std::string &filename);
void file_vtk_pvd_reformat(data &d);

#endif /* file_hpp */
