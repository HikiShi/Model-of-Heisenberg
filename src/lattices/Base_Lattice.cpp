//
// Created by hiki on 3/21/19.
//

#include "Base_Lattice.h"

// TODO: Add initializer for different lattices types
//Currently, support only ferromagnetic lattices

Base_Lattice::Base_Lattice(int const& size_of_x, int const& size_of_y, int const& size_of_z, spin_type const& J) : size_of_x{size_of_x}, size_of_y{size_of_y}, size_of_z{size_of_z}, J{J}
{
    random_z = std::uniform_real_distribution<>(-1 ,1);
    random_teta = std::uniform_real_distribution<>(0, 2 * 3.14159265358979323846);

    if(size_of_x <= 0 || size_of_y <= 0 || size_of_z <= 0)
    {
        throw 1;
    }
    lattice = new Spin** [size_of_x];
    for( int index_x = 0; index_x < size_of_x; index_x++)
    {
        lattice[index_x] = new Spin* [size_of_y];
        for( int index_y = 0; index_y < size_of_y; index_y++)
        {
            lattice[index_x][index_y] = new Spin [size_of_z];
        }
    }
}

spin_type Base_Lattice::get_energy_of_two_spin(Spin spin1, Spin spin2)
{
    return -J  * ((spin1.x * spin2.x) + (spin1.y * spin2.y) + (spin1.z * spin2.z));
}


spin_type Base_Lattice::get_energy_spin_and_neighbors(int const& index_x, int const& index_y, int const& index_z)
{
    spin_type energy = 0;

    index_x == 0
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[size_of_x - 1][index_y][index_z])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x - 1][index_y][index_z]);

    index_x == size_of_x - 1
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[0][index_y][index_z])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x + 1][index_y][index_z]);

    index_y == 0
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][size_of_y - 1][index_z])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y - 1][index_z]);

    index_y == size_of_y - 1
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][0][index_z])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y + 1][index_z]);

    index_z == 0
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y][size_of_z - 1])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y][index_z - 1]);

    index_z == size_of_z - 1
        ? energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y][0])
        : energy += get_energy_of_two_spin(lattice[index_x][index_y][index_z], lattice[index_x][index_y][index_z + 1]);

    return energy;
}

out_data Base_Lattice::get_lattice_date()
{
    out_data data;
    spin_type energy = 0.0;
    spin_type x_components = 0.0;
    spin_type y_components = 0.0;
    spin_type z_components = 0.0;

    for (int i = 0; i < size_of_x; i++)
    {
        for (int j = 0; j < size_of_y; j++)
        {
            for (int k = 0; k < size_of_z; k++)
            {
                energy += get_energy_spin_and_neighbors(i, j, k);

                x_components += lattice[i][j][k].x;
                y_components += lattice[i][j][k].y;
                z_components += lattice[i][j][k].z;
            }
        }
    }

    data.e += energy;
    data.m += fabs(static_cast<spin_type>(1.0 / (size_of_x * size_of_y * size_of_z)) * pow((pow(x_components, 2.0) + pow(y_components, 2.0) + pow(z_components, 2.0)), 0.5));

    return data;
}

 void Base_Lattice::set_new_configuration(int const& index_x, int const& index_y, int const& index_z)
 {
     Spin random_spin;

     random_spin.z = random_z(generator);
     const spin_type teta = random_teta(generator);
     const spin_type temp = sqrt(1.0 - (random_spin.z * random_spin.z));
     random_spin.x = temp * cos(teta);
     random_spin.y = temp * sin(teta);

     lattice[index_x][index_y][index_z] = random_spin;
 }

Spin Base_Lattice::get_spin(int const& index_x, int const& index_y, int const& index_z)
{
    return lattice[index_x][index_y][index_z];
}

void Base_Lattice::set_spin_on_position(Spin const& spin, int const& index_x, int const& index_y, int const& index_z)
{
    lattice[index_x][index_y][index_z] = spin;
}

Base_Lattice::~Base_Lattice()
{
    for(int index_x = 0; index_x <size_of_x; index_x++)
    {
        for(int index_y = 0; index_y < size_of_y; index_y++)
        {
            delete[] lattice[index_x][index_y];
        }
        delete[] lattice[index_x];
    }
    delete[] lattice;
};

//
//template <typename spin_type>
//Lattice_Proxy_2<spin_type>::Lattice_Proxy_2(Spin *Lattice_row)  : Lattice_row{Lattice_row} {}
//
//template <typename spin_type>
//Spin& Lattice_Proxy_2<spin_type>::operator[](int index)
//{
//    return Lattice_row[index];
//}
//
//template <typename  spin_type>
//Lattice_Proxy<spin_type>::Lattice_Proxy(Spin **Lattice_Plane) : Lattice_Plane{Lattice_Plane} {}
//
//template <typename  spin_type>
//Lattice_Proxy_2<spin_type> Lattice_Proxy<spin_type>::operator[](int index)
//{
//    return Lattice_Proxy_2(Lattice_Plane[index]);
//}
//
//template <typename spin_type>
//Lattice_Proxy<spin_type>  Lattice<spin_type>::operator[](int index)
//{
//    return Lattice_Proxy(_Lattice[index]);
//}