//
//  update.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "update.hpp"

void update(data &d)
{
    d.fluid.vel_o = d.fluid.vel;
    d.fluid.vel = d.fluid.vel_n;
    d.fluid.pos_o = d.fluid.pos;
    d.fluid.pos = d.fluid.pos_n;
    d.fluid.rho_o = d.fluid.rho;
    d.fluid.rho = d.fluid.rho_n;
    d.time.dt_o = d.time.dt;
    d.time.t += d.time.dt;                                        // set new time
}
