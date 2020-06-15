//
//  main.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include <iostream>
#include <chrono>
#include <cmath>
#include "data.h"
#include "init.hpp"
#include "setup.hpp"
#include "solve.hpp"
#include "output.hpp"
#include "update.hpp"
#include "utils.hpp"

int main(int argc, const char * argv[])
{
    data d;                                                 // define data structure
    
    std::cout << std::endl << argv[0] << " --- created by Yeaw Chu Lee, 2015" << std::endl;
    
    init(d, argc, argv);                                    // initialise solver
    
    timer<std::chrono::high_resolution_clock> stopwatch;    // define stopwatch
    stopwatch.start();                                      // start stopwatch

    long itr = 0;                                           // iteration number
    while ((d.time.t < d.time.t_end) || near_equal(d.time.t, d.time.t_end))
    {
        setup(d);                                           // setup and condition variables
        solve(d);                                           // solve variables
        output(d);                                          // output variables
        update(d);                                          // update variables
        ++itr;                                              // increase iteration count
    }
        
    stopwatch.stop();                                       // stop stopwatch
    std::cout << std::endl;
    std::cout << "Elapse time: " << stopwatch << std::endl;
    std::cout << "Total iterations: " << itr << std::endl;
    
    return 0;                                               // return error code 0
}
