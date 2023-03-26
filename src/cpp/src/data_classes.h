#ifndef DATA_CLASSES_H
#define DATA_CLASSES_H

#include <string>

class Fluid
{

public:
    std::string Name;
    double CriticalPressure;
    double CriticalTemperature;
    double CriticalCompressibility;
    double AccentricFactor;
    double LennardJonnesDiameter;
    double LennardJonnesEnergy;
    ~Fluid() {}
    Fluid() {}
    Fluid(std::string name, double critical_pressure, double critical_temperature, double accentric_factor, double critical_compressibility = 0.0, double lj_diameter = 0.0, double lj_energy = 0.0)
    {
        this->Name = name;
        this->CriticalPressure = critical_pressure;
        this->CriticalTemperature = critical_temperature;
        this->AccentricFactor = accentric_factor;
        this->CriticalCompressibility = critical_compressibility;
        this->LennardJonnesDiameter = lj_diameter;
        this->LennardJonnesEnergy = lj_energy;
    }
};

class Adsorbent
{

public:
    std::string Name;
    double SolidDiameter;      // sigma_ss, Angstrom
    double SolidAtomicDensity; // rho_atoms, Angstrom
    Adsorbent() {}
    Adsorbent(std::string name, double solid_diameter, double solid_atomic_density)
    {
        this->Name = name;
        this->SolidDiameter = solid_diameter;
        this->SolidAtomicDensity = solid_atomic_density;
    }

    ~Adsorbent() {}
};

#endif