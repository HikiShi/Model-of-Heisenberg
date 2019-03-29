cmake_minimum_required(VERSION 3.10.2)
project(Model_of_Heisenberg)

set(CMAKE_CXX_STANDARD 17)

add_executable( Model_of_Heisenberg src/model.cpp src/components_of_Hamiltonian.h src/spin.h src/configuration_of_system.h)