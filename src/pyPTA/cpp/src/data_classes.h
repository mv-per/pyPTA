// some_header_file.h
#ifndef DATA_CLASSES_H
#define DATA_CLASSES_H

struct Fluid
{

    Fluid(){};
    ~Fluid(){};

    double critical_pressure;
    double critical_temperature;
    double critical_compressibility;
    double acentric_factor;
    double lj_diameter;
    double lj_energy;
};

struct Adsorbent
{

    Adsorbent() {}
    Adsorbent(double sigma_ss_, double rho_atoms_)
    {
        sigma_ss = sigma_ss_;
        rho_atoms = rho_atoms_;
    }
    ~Adsorbent() {}

    double sigma_ss, rho_atoms;
};

#endif