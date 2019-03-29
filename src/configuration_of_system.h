//
// Created by lotos on 29.03.2019.
//

#ifndef MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H
#define MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H


constexpr int L = 16;
constexpr int N = 3;
constexpr int count_of_layers = 2;
constexpr double J = 1.0; //integral of interaction
constexpr double J2 = - 0.3 * J;
constexpr int time_of_relax = 1000;
constexpr int time_of_observation = 10000;
constexpr int temperature_steps_count = 30;
constexpr int static_configuration_count = 1;
constexpr double delta = 0.68; //const of anisotropy
constexpr double temperature_start = 0.8;
constexpr double size_of_temperature_step = 0.05;


#endif //MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H
