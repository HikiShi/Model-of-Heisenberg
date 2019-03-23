cmake_minimum_required(VERSION 3.10.2)
project(Model_of_Heisenberg)

set(CMAKE_CXX_STANDARD 17)

add_executable(
        Model_of_Heisenberg
        src/main.cpp
        src/lattice_builder/Lattice_builder.h
        src/lattices/Base_Lattice.h
        src/lattices/Spin.h
        src/algorithms/metropolis_algorithm.h src/output/output_data.h)