#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include "components_of_Hamiltonian.h"
#include <unordered_set>

using namespace std;
using namespace chrono;

mt19937_64 generator;
uniform_real_distribution<double> z_rand(-1, 1);
uniform_real_distribution<double> teta_rand(0, 2 * M_PI);

unordered_set<int> rows_with_upper_interaction;
unordered_set<int> rows_with_lower_interaction;


inline spin get_random_spin() noexcept
{
    spin random_spin;

    random_spin.z = z_rand(generator);
    const double teta = teta_rand(generator);
    const double temp = sqrt(1 - random_spin.z * random_spin.z);
    random_spin.x = temp * cos(teta);
    random_spin.y = temp * sin(teta);

    return random_spin;
}

struct out_data
{
    double e = 0;
    double m = 0;
    double e2 = 0;
    double m2 = 0;
};

inline void update_out_data(spin const *const *const *const sp, out_data &data) noexcept
{
    double x_components = 0.0;
    double y_components = 0.0;
    double z_components = 0.0;
    double energy = 0.0;

    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N * count_of_layers; k++)
            {
                energy += get_energy_of_spin_and_neighbours(sp, rows_with_upper_interaction, rows_with_lower_interaction,  i, j, k);

                x_components += sp[i][j][k].x;
                y_components += sp[i][j][k].y;
                z_components += sp[i][j][k].z;
            }
        }
    }
    data.e += energy;
    data.e2 += pow(energy, 2);
    data.m += fabs(static_cast<double>(1.0 / (N * count_of_layers * L * L)) * pow((pow(x_components, 2.0) + pow(y_components, 2.0) + pow(z_components, 2.0)), 0.5));
    data.m2 += pow(static_cast<double>(1.0 / (N * count_of_layers * L * L)) * pow((pow(x_components, 2.0) + pow(y_components, 2.0) + pow(z_components, 2.0)), 0.5), 2);
}

inline  void normalize_out_data(out_data &data, int const &count_of_steps)
{
    data.e  /= static_cast<double>(count_of_steps);
    data.m  /= static_cast<double>(count_of_steps);
    data.e2 /= static_cast<double>(count_of_steps);
    data.m2 /= static_cast<double>(count_of_steps);
}


int main()
{

    auto time_of_start = high_resolution_clock::now();
    for(int i = 1; i < count_of_layers; i++)
    {
        rows_with_upper_interaction.insert(i * N);
        rows_with_lower_interaction.insert((i * N) - 1);
    }

    cout << rows_with_lower_interaction.size() << endl;
    cout << rows_with_upper_interaction.size() << endl;

    out_data data_for_file[temperature_steps_count];

    for (int static_configuration = 0; static_configuration < static_configuration_count; static_configuration++)
    {
        spin ***lattice = new spin **[L];
        for (size_t i = 0; i < L; i++)
        {
            lattice[i] = new spin *[L];
            for (size_t j = 0; j < L; j++)
            {
                lattice[i][j] = new spin[N * count_of_layers];
            }
        }

        cout << static_configuration << endl;
        uniform_real_distribution<double> get_r(0, 1);
        uniform_int_distribution<int> get_i_or_j(0, L - 1);
        uniform_int_distribution<int> get_k(0, N * count_of_layers - 1);

        for (int temperature_step = 0; temperature_step < temperature_steps_count; temperature_step++)
        {
            const double temperature = temperature_start + (temperature_step * size_of_temperature_step);
            int time_of_relaxation = temperature_step == 0 ? 1000 : time_of_relax;
            int time_of_simulation = time_of_relaxation + time_of_observation;

            int count_of_iterations = 0;
            for (int iteration = 0; iteration < time_of_simulation; iteration++)
            {
              
                if (iteration >= time_of_relaxation)
                {
                    update_out_data(lattice, data_for_file[temperature_step]);
                }

                for (int elem_step = 0; elem_step < N * count_of_layers * L * L; ++elem_step)
                {
                    int i = get_i_or_j(generator);
                    int j = get_i_or_j(generator);
                    int k = get_k(generator);

                    const spin selected_spin = lattice[i][j][k];
                    const double E1 = get_energy_of_spin_and_neighbours(lattice, rows_with_upper_interaction, rows_with_lower_interaction,  i, j, k);
                    lattice[i][j][k] = get_random_spin();

                    const double E2 = get_energy_of_spin_and_neighbours(lattice, rows_with_upper_interaction, rows_with_lower_interaction,  i, j, k);

                    const double delta_E =  E2 - E1;


                    const double W = exp(- delta_E / ((temperature)));

                    if (delta_E > 0)
                    {
                        const double rand_value = get_r(generator);
                        if (rand_value > W)
                        {
                            lattice[i][j][k] = selected_spin;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                delete[] lattice[i][j];
            }
            delete[] lattice[i];
        }
        delete[] lattice;
    }

    for (int temperature_step = 0; temperature_step < temperature_steps_count; temperature_step++)
    {
        normalize_out_data(data_for_file[temperature_step], time_of_observation * static_configuration_count);
    }

    std::ofstream output_file("output_data.dat", std::ios_base::trunc);

    for (int temperature_step = 0; temperature_step < temperature_steps_count; temperature_step++)
    {
        const double temperature = temperature_start + (temperature_step * size_of_temperature_step);
        const double m = data_for_file[temperature_step].m;
        const double e = data_for_file[temperature_step].e;

        const double m_sqr = data_for_file[temperature_step].m2;
        const double chi = L * L * N * count_of_layers * (m_sqr - m * m) / temperature;

        const double e_sqr = data_for_file[temperature_step].e2;
        const double C = L * L * N * count_of_layers * (e_sqr - e * e) / (temperature * temperature);

        output_file << std::fixed << std::setprecision(15) << temperature << "\t" << m << "\t" << e << "\t" << chi << "\t" << C << std::endl;
    }

    output_file.close();

    auto time_of_end = high_resolution_clock::now();

    cout << "I SPEND " << duration_cast<seconds>(time_of_end - time_of_start).count() << " SECONDs OF MY LIFE" << endl;
}
