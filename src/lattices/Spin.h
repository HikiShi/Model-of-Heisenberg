//
// Created by hiki on 3/21/19.
//
#ifndef MODEL_OF_HEISENBERG_SPIN_H
#define MODEL_OF_HEISENBERG_SPIN_H

using spin_type = double;

struct Spin {
    spin_type x = 0.0;
    spin_type y = 0.0;
    spin_type z = 1.0;
};

#endif //MODEL_OF_HEISENBERG_SPIN_H