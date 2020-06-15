//
//  init.hpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef init_hpp
#define init_hpp

#include <iostream>
#include "data.h"
#include "file.hpp"
#include "kernel.hpp"

void init(data &d, int argc, const char * argv[]);
void init_config(data &d);
void init_array(data &d);
void init_kernel(data &d);

#endif /* init_hpp */
