#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include <string>
#include "lattices/Base_Lattice.h"
#include "algorithms/metropolis_algorithm.h"
#include <chrono>

int main()
{
   auto t_1 = std::chrono::high_resolution_clock::now();
   auto lattice = Base_Lattice(16, 16, 3, 1.0);
   metropolis_algorithm(lattice, 10, 10000, 1000, 0.7, 50, 0.05);
   auto t_2 = std::chrono::high_resolution_clock::now();

   std::cout << "spend seconds" << std::chrono::duration_cast<std::chrono::seconds>(t_2 - t_1).count() << std::endl;
}
