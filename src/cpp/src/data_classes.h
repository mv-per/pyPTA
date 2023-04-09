#ifndef DATA_CLASSES_H
#define DATA_CLASSES_H

#include <string>
#include <optional>

class Fluid
{

public:
    std::string Name;
    double CriticalPressure;
    double CriticalTemperature;
    std::optional<double> CriticalCompressibility;
    double AccentricFactor;
    std::optional<double> LennardJonnesDiameter;
    std::optional<double> MolecularWeight;
    ~Fluid() {}
    Fluid() {}
    Fluid(std::string Name,
          double CriticalPressure,
          double CriticalTemperature,
          double AccentricFactor,
          double CriticalCompressibility = 0.0,
          double LennardJonnesDiameter = 0.0,
          double MolecularWeight = 0.0)
    {
        this->Name = Name;
        this->CriticalPressure = CriticalPressure;
        this->CriticalTemperature = CriticalTemperature;
        this->AccentricFactor = AccentricFactor;
        this->CriticalCompressibility = CriticalCompressibility;
        this->LennardJonnesDiameter = LennardJonnesDiameter;
        this->MolecularWeight = MolecularWeight;
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