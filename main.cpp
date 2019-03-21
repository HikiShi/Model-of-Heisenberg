#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <cassert>
#include <unordered_set>

using namespace std;
using namespace chrono;

constexpr int L = 16;
constexpr int N = 3;
constexpr double J = 1.0; //integral of interaction
constexpr int time_of_relax = 1000;
constexpr int time_of_observation = 10000;
constexpr int temperature_steps_count = 200;
constexpr int static_configuration_count = 10;
constexpr double delta = 0.66; //const of anisotropy

mt19937_64 generator;
uniform_real_distribution<double> z_rand(-1, 1);
uniform_real_distribution<double> teta_rand(0, 2 * M_PI);

double get_random_z()
{

    return z_rand(generator);
}

double get_random_teta()
{

    return teta_rand(generator);
}

struct spin
{
    double x{0.0};
    double y{0.0};
    double z{1.0};
};

spin get_random_spin()
{
    spin random_spin;

    random_spin.z = get_random_z();
    const double teta = get_random_teta();
    const double temp = sqrt(1 - random_spin.z * random_spin.z);
    random_spin.x = temp * cos(teta);
    random_spin.y = temp * sin(teta);

    return random_spin;
}

inline double get_energy_of_two_spin(spin const &spin1, spin const &spin2)
{
    double temp = ((spin1.x * spin2.x) + (spin1.y * spin2.y));
    temp = (1.0 - delta) * temp;
    temp = temp + (spin1.z * spin2.z);
    return -J * temp;
}

double get_energy_of_spin_and_neighbours(spin const *const *const *const sp, int const &i, int const &j, int const &k)
{
    double energy = 0;
    if (i == 0 || j == 0 || k == 0 || i == L - 1 || j == L - 1 || k == N - 1)
    {
        //make marginal
        i == 0
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[L - 1][j][k])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i - 1][j][k]);

        i == L - 1
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[0][j][k])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i + 1][j][k]);

        j == 0
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][L - 1][k])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j - 1][k]);

        j == L - 1
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][0][k])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j + 1][k]);

        k == 0
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][N - 1])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);

        k == N - 1
            ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][0])
            : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);

        return energy;
    }

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i + 1][j][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i - 1][j][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j + 1][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j - 1][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);

    return energy;
}

struct out_data
{
    double e = 0;
    double m = 0;
    double e2 = 0;
    double m2 = 0;
};

struct selected_spin
{
    int i = 0;
    int j = 0;
    int k = 0;
};

void update_out_data(spin const *const *const *const sp, out_data &data)
{
    double x_components = 0.0;
    double y_components = 0.0;
    double z_components = 0.0;
    double energy = 0.0;

    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N; k++)
            {
                energy += get_energy_of_spin_and_neighbours(sp, i, j, k);

                x_components += sp[i][j][k].x;
                y_components += sp[i][j][k].y;
                z_components += sp[i][j][k].z;
            }
        }
    }
    data.e += energy;
    data.e2 += pow(energy, 2);
    data.m += fabs(static_cast<double>(1.0 / (N * L * L)) * pow((pow(x_components, 2.0) + pow(y_components, 2.0) + pow(z_components, 2.0)), 0.5));
    data.m2 += pow(static_cast<double>(1.0 / (N * L * L)) * pow((pow(x_components, 2.0) + pow(y_components, 2.0) + pow(z_components, 2.0)), 0.5), 2);
}

void normalize_out_data(out_data &data, int const &count_of_steps)
{
    data.e /= static_cast<double>(count_of_steps);
    data.m /= static_cast<double>(count_of_steps);
    data.e2 /= static_cast<double>(count_of_steps);
    data.m2 /= static_cast<double>(count_of_steps);
}

mt19937_64 random_index;
uniform_real_distribution<double> get_r(0,
                                        1);
uniform_int_distribution<int> get_i_or_j(0, L - 1);
uniform_int_distribution<int> get_k(0, N - 1);

// void f1(spin ***lattice, out_data *data_for_file, int &count_of_interations, int &time_of_simulation, int &temperature_step, int &time_of_relaxation, mutex &mutex_count, mutex &mutex_out_data) //, )
// {

//     

// }

int main()
{

    auto time_of_start = high_resolution_clock::now();

    out_data data_for_file[temperature_steps_count];

    for (int static_configuration = 0; static_configuration < static_configuration_count; static_configuration++)
    {
        spin ***lattice = new spin **[L];
        for (size_t i = 0; i < L; i++)
        {
            lattice[i] = new spin *[L];
            for (size_t j = 0; j < L; j++)
            {
                lattice[i][j] = new spin[N];
            }
        }

        cout << static_configuration << endl;
        mt19937_64 random_index;
        uniform_real_distribution<double> get_r(0, 1);
        uniform_int_distribution<int> get_i_or_j(0, L - 1);
        uniform_int_distribution<int> get_k(0, N - 1);

        for (int temperature_step = 0; temperature_step < temperature_steps_count; temperature_step++)
        {
            const double temperature = 2.08 + (temperature_step * 0.02);
            int time_of_relaxation = temperature_step == 0 ? 10000 : time_of_relax;
            int time_of_simulation = time_of_relaxation + time_of_observation;

            int count_of_iterations = 0;
            for (int iteration = 0; iteration < time_of_simulation; iteration++)
            {
              
                if (iteration >= time_of_relaxation)
                {
                    update_out_data(lattice, data_for_file[temperature_step]);
                }

                for (int elem_step = 0; elem_step < N * L * L; ++elem_step)
                {
                    int i = get_i_or_j(random_index);
                    int j = get_i_or_j(random_index);
                    int k = get_k(random_index);

                    const spin selected_spin = lattice[i][j][k];
                    const double E1 = get_energy_of_spin_and_neighbours(lattice, i, j, k);
                    lattice[i][j][k] = get_random_spin();

                    const double E2 = get_energy_of_spin_and_neighbours(lattice, i, j, k);

                    const double delta_E =  E2 - E1;


                    const double W = exp(- delta_E / ((temperature)));

                    if (delta_E > 0)
                    {
                        //cout << iteration << " " << elem_step << "  " << delta_E << endl;
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
        const double temperature = 2.08 + (temperature_step * 0.02);
        const double m = data_for_file[temperature_step].m;
        const double e = data_for_file[temperature_step].e;

        const double m_sqr = data_for_file[temperature_step].m2;
        const double chi = L * L * N * (m_sqr - m * m) / temperature;

        const double e_sqr = data_for_file[temperature_step].e2;
        const double C = L * L * N * (e_sqr - e * e) / (temperature * temperature);

        output_file << std::fixed << std::setprecision(15) << temperature << "\t" << m << "\t" << e << "\t" << chi << "\t" << C << std::endl;
    }

    output_file.close();

    auto time_of_end = high_resolution_clock::now();

    cout << "I SPEND " << duration_cast<seconds>(time_of_end - time_of_start).count() << " SECONDs OF MY LIFE" << endl;
}
