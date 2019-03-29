//
// Created by lotos on 29.03.2019.
//

#ifndef MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H
#define MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H

#include "spin.h"
#include "configuration_of_system.h"
#include <unordered_set>


inline double get_energy_of_two_spin(spin const &spin1, spin const &spin2) noexcept
{
    return -J *  (((1.0 - delta) * ((spin1.x * spin2.x) + (spin1.y * spin2.y)) ) +   (  spin1.z * spin2.z)  );

}

inline double get_energy_of_two_spin_whith_rkky_interaction(spin const &spin1, spin const &spin2) noexcept
{
    return  -J *  (((1.0 - delta) * ((spin1.x * spin2.x) + (spin1.y * spin2.y))) + (spin1.z * spin2.z)) - J2 * (((1.0 - delta) * ((spin1.x * spin2.x) + (spin1.y * spin2.y))) +   (spin1.z * spin2.z));

}



inline double get_energy_of_spin_and_neighbours_with_uppder_interaction(spin const *const *const *const sp, int const &i, int const &j, int const &k) noexcept
{
    double energy = 0;

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

    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i][j][k - 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);

    return energy;
}

inline double get_energy_of_spin_and_neighbours_with_lower_interaction(spin const *const *const *const sp, int const &i, int const &j, int const &k) noexcept
{
    double energy = 0;

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

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i][j][k + 1]);

    return energy;
}


inline double get_energy_of_spin_and_neighbours_without_rrky_interaction(spin const *const *const *const sp, int const &i, int const &j, int const &k) noexcept
{
    double energy = 0;


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

    if(k != 0)
    {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);
    }

    if(k != N * count_of_layers - 1)
    {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);
    }

    return energy;
}

inline double get_energy_of_spin_and_neighbours (spin const *const *const *const sp, std::unordered_set<int> const& rows_with_upper_interaction,
        std::unordered_set<int> const& rows_with_lower_interaction,  int const &i, int const &j, int const &k) noexcept
{
    if(rows_with_upper_interaction.find(k) != rows_with_upper_interaction.end())
    {
        return get_energy_of_spin_and_neighbours_with_uppder_interaction(sp, i, j, k);
    }

    if(rows_with_lower_interaction.find(k) != rows_with_lower_interaction.end())
    {
        return get_energy_of_spin_and_neighbours_with_lower_interaction(sp, i, j, k);
    }

    return get_energy_of_spin_and_neighbours_without_rrky_interaction(sp, i, j, k);

}




#endif //MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H
