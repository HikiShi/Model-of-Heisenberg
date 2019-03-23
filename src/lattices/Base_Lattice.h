//
// Created by hiki on 3/21/19.
//

#ifndef MODEL_OF_HEISENBERG_LATTICE_H
#define MODEL_OF_HEISENBERG_LATTICE_H

#include "Spin.h"
#include <random>
#include <iostream>
#include "../output/output_data.h"


class Base_Lattice {
private:
    std::mt19937 generator;
    std::uniform_real_distribution<> random_z;
    std::uniform_real_distribution<> random_teta;

public:
    Spin*** lattice;
    int size_of_x;
    int size_of_y;
    int size_of_z;
    spin_type J;

    Base_Lattice(int const&, int const&, int const&, spin_type const&);
    virtual spin_type get_energy_of_two_spin(Spin, Spin);
    virtual spin_type get_energy_spin_and_neighbors(int const&, int const&, int const&);
    virtual out_data get_lattice_date();
    virtual Spin get_spin(int const&, int const&, int const&);
    virtual void set_spin_on_position(Spin const&, int const&, int const&, int const&);
    virtual void set_new_configuration(int const&, int const&, int const&);
    virtual ~Base_Lattice();
};





//class Lattice {
//private:
//    Spin*** lattices;
//    double J;
//    double delta;
//public:
//    int size_of_x;
//    int size_of_y;
//    int size_of_z;
//
//    Lattice();
//    ~Lattice();
//};





#include "Base_Lattice.cpp"
#endif //MODEL_OF_HEISENBERG_LATTICE_H
