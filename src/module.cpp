

#include "cpp/src/equations_of_state/eos.h"
#include "cpp/src/equations_of_state/eos_helper.h"
#include "cpp/src/equations_of_state/pr77.h"
#include "cpp/src/optimization_algorithms/brent.h"
#include "cpp/src/data_classes.h"
#include "cpp/src/adsorption_potentials.h"
#include "cpp/src/helpers.h"
#include "cpp/src/pta_helper.h"
#include "cpp/src/pta_solver.h"
#include "cpp/src/pta_pure.h"
#include "cpp/src/pta_mixture.h"

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pyPTA, m)
{
    m.doc() = "C/C++ methods to run the PTA model";

    py::class_<Adsorbent>(m, "Adsorbent")
        .def(py::init<std::string &, double &, double &>())
        .def_readwrite("Name", &Adsorbent::Name)
        .def_readwrite("SolidDiameter", &Adsorbent::SolidDiameter)
        .def_readwrite("SolidAtomicDensity", &Adsorbent::SolidAtomicDensity);

    py::class_<Fluid>(m, "Fluid")
        .def(py::init<>())
        .def_readwrite("CriticalPressure", &Fluid::CriticalPressure)
        .def_readwrite("CriticalTemperature", &Fluid::CriticalTemperature)
        .def_readwrite("CriticalCompressibility", &Fluid::CriticalCompressibility)
        .def_readwrite("AccentricFactor", &Fluid::AccentricFactor)
        .def_readwrite("LennardJonnesDiameter", &Fluid::LennardJonnesDiameter)
        .def_readwrite("LennardJonnesEnergy", &Fluid::LennardJonnesEnergy);

    py::class_<PurePTA>(m, "PurePTA")
        .def(py::init<std::string &, std::string &, std::string &, std::size_t &>(),
             py::arg("adsorption_potential"),
             py::arg("equation_of_state"),
             py::arg("isotherm_type"),
             py::arg("number_of_layers"))
        .def_readwrite("Potential", &PurePTA::Potential)
        .def_readwrite("EquationOfState", &PurePTA::EquationOfState)
        .def_readwrite("IsothermType", &PurePTA::IsothermType)
        .def_readwrite("NumberOfLayers", &PurePTA::NumberOfLayers)
        .def("GetLoading", &PurePTA::GetLoading,
             "Get single-component loading",
             py::arg("pressure"),
             py::arg("temperature"),
             py::arg("adsorption_potential_parameters"),
             py::arg("fluid"))
        .def("GetDeviationRange", &PurePTA::GetDeviationRange,
             "Get absolute difference from experimental and calculated loadings",
             py::arg("deviation_type"),
             py::arg("experimental_loading"),
             py::arg("experimental_pressure"),
             py::arg("temperature"),
             py::arg("adsorption_potential_parameters"),
             py::arg("fluid"))
        .def("GetLoadings", &PurePTA::GetLoadings,
             "Get single-component loading for multiple pressure datapoints",
             py::arg("experimental_pressure"),
             py::arg("temperature"),
             py::arg("adsorption_potential_parameters"),
             py::arg("fluid"));

    // py::class_<MixturePTA>(m, "MixturePTA")
    //     .def(py::init<std::string &, std::string &, std::string &, std::size_t &>(),
    //          py::arg("adsorption_potential"),
    //          py::arg("equation_of_state"),
    //          py::arg("isotherm_type"),
    //          py::arg("number_of_layers"))
    //     .def_readwrite("adsorption_potential", &MixturePTA::potential)
    //     .def_readwrite("equation_of_state", &MixturePTA::equation_of_state)
    //     .def_readwrite("isotherm_type", &MixturePTA::isotherm_type)
    //     .def_readwrite("number_of_layers", &MixturePTA::num_of_layers)
    //     .def("GetLoading", &MixturePTA::GetLoading,
    //          "Get multiple-component loadings",
    //          py::arg("composition"),
    //          py::arg("pressure"),
    //          py::arg("temperature"),
    //          py::arg("adsorption_potential_parameters"),
    //          py::arg("fluid_parameters"));
}
