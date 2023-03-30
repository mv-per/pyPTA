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
    Fluid(std::string Name,
          double CriticalPressure,
          double CriticalTemperature,
          double AccentricFactor,
          double CriticalCompressibility = 0.0,
          double LennardJonnesDiameter = 0.0,
          double LennardJonnesEnergy = 0.0)
    {
        this->Name = Name;
        this->CriticalPressure = CriticalPressure;
        this->CriticalTemperature = CriticalTemperature;
        this->AccentricFactor = AccentricFactor;
        this->CriticalCompressibility = CriticalCompressibility;
        this->LennardJonnesDiameter = LennardJonnesDiameter;
        this->LennardJonnesEnergy = LennardJonnesEnergy;
    }
};

class Adsorbent
{

public:
    std::string Name;
    double SolidDiameter;      // sigma_ss, Angstrom
    double SolidAtomicDensity; // rho_atoms, Angstrom
    Adsorbent() {}
    Adsorbent(std::string Name, double SolidDiameter, double SolidAtomicDensity)
    {
        this->Name = Name;
        this->SolidDiameter = SolidDiameter;
        this->SolidAtomicDensity = SolidAtomicDensity;
    }

    ~Adsorbent() {}
};

const std::string DRA_POTENTIAL = "DRA";
const std::string LEE_POTENTIAL = "LEE";
const std::string STEELE_POTENTIAL = "STEELE";

#endif