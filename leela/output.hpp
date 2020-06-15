//
//  output.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef output_hpp
#define output_hpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include "data.h"
#include "file.hpp"
#include "utils.hpp"

void output(data &d);
void output_file(data &d);
void output_screen(data &d);
void output_screen_progress(data &d);
void output_screen_progress_info(data &d);
void output_screen_progress_bar(data &d);
void output_screen_progress_minimal(data &d);
void output_screen_progress_debug(data &d);

#endif /* output_hpp */
