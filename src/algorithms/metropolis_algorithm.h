//
// Created by lotos on 26.03.2019.
//

#ifndef MODEL_OF_HEISENBERG_METROPOLIS_ALGORITHM_H
#define MODEL_OF_HEISENBERG_METROPOLIS_ALGORITHM_H

#include <random>
#include "../lattices/Base_Lattice.h"
#include "../lattices/Spin.h"
#include "../output/output_data.h"
#include <fstream>
#include <iostream>
#include <iomanip>



inline void update_out_data(out_data data_from_lattice, out_data& data)
{
    data.e += data_from_lattice.e;
    data.m += data_from_lattice.m;
}


void metropolis_algorithm(Base_Lattice& lattice,
                          int const& count_of_static_configuration,
                          int const& time_of_observation,
                          int const& time_of_relaxation,
                          double const& temperature_start,
                          int const& temperature_step_count,
                          double temperature_step_size)
{
    out_data data_for_file[temperature_step_count];
    for(int configuration = 0; configuration < count_of_static_configuration; configuration++)
    {
        std::cout << configuration << std::endl;

        std::mt19937_64 random_index(configuration);
        std::uniform_real_distribution<double> get_r(0, 1);
        std::uniform_int_distribution<int> get_i_or_j(0, lattice.size_of_x - 1);
        std::uniform_int_distribution<int> get_k(0, lattice.size_of_z - 1);


        for (int temperature_step = 0; temperature_step < temperature_step_count; temperature_step++) {
            const double temperature = temperature_start + (temperature_step * temperature_step_size);
            int time_of_simulation = time_of_relaxation + time_of_observation;

            for (int iteration = 0; iteration < time_of_simulation; iteration++) {
                if (iteration >= time_of_relaxation) {
                    update_out_data(lattice.get_lattice_date(), data_for_file[temperature_step]);
                }

                for (int elem_step = 0;
                     elem_step < lattice.size_of_z * lattice.size_of_x * lattice.size_of_y; ++elem_step) {
                    int i = get_i_or_j(random_index);
                    int j = get_i_or_j(random_index);
                    int k = get_k(random_index);

                    const Spin selected_spin = lattice.get_spin(i, j, k);
                    const double E1 = lattice.get_energy_spin_and_neighbors(i, j, k);
                    lattice.set_new_configuration(i, j, k);

                    const double E2 = lattice.get_energy_spin_and_neighbors(i, j, k);

                    const double delta_E = E2 - E1;


                    const double W = exp(-delta_E / ((temperature)));

                    if (delta_E > 0) {
                        const double rand_value = get_r(random_index);
                        if (rand_value > W) {
                            lattice.set_spin_on_position(selected_spin, i, j, k);
                        }
                    }
                }
            }
        }
    }


    for (int temperature_step = 0; temperature_step < temperature_step_count; temperature_step++)
    {
        data_for_file[temperature_step].e /= static_cast<spin_type>(time_of_observation * count_of_static_configuration);
        data_for_file[temperature_step].m /= static_cast<spin_type>(time_of_observation * count_of_static_configuration);
    }

    std::ofstream output_file("output_data.dat", std::ios_base::trunc);

    for (int temperature_step = 0; temperature_step < temperature_step_count; temperature_step++)
    {
        const double temperature = temperature_start + (temperature_step * temperature_step_size);
        const double e = data_for_file[temperature_step].e;
        const double m = data_for_file[temperature_step].m;
        output_file << std::fixed << std::setprecision(15) << temperature << "\t" << m << "\t" << e << std::endl;
    }

    output_file.close();

}




#endif //MODEL_OF_HEISENBERG_METROPOLIS_ALGORITHM_H
