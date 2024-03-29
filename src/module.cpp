

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "cpp/src/pta_pure.h"
#include "cpp/src/pta_mixture.h"

PYBIND11_MODULE(pyPTA, m)
{
     m.doc() = "C/C++ methods to run the PTA model";

     py::class_<Adsorbent>(m, "Adsorbent")
         .def(py::init<>())
         .def_readwrite("name", &Adsorbent::Name)
         .def_readwrite("diameter", &Adsorbent::SolidDiameter)
         .def_readwrite("atomic_density", &Adsorbent::SolidAtomicDensity)
         .def(py::init([](std::string name, double solid_diameter, double solid_atomic_density)
                       { return new Adsorbent(name, solid_diameter, solid_atomic_density); }),
              py::arg("name"),
              py::arg("diameter"),
              py::arg("atomic_density"));

     py::class_<Fluid>(m, "Fluid")
         .def_readwrite("name", &Fluid::Name)
         .def_readwrite("critical_pressure", &Fluid::CriticalPressure)
         .def_readwrite("critical_temperature", &Fluid::CriticalTemperature)
         .def_readwrite("accentric_factor", &Fluid::AccentricFactor)
         .def_readwrite("critical_compressibility", &Fluid::CriticalCompressibility)
         .def_readwrite("lennard_jonnes_diameter", &Fluid::LennardJonnesDiameter)
         .def(py::init([](std::string name,
                          double critical_pressure,
                          double critical_temperature,
                          double accentric_factor,
                          double critical_compressibility,
                          double lennard_jonnes_diameter)
                       { return new Fluid(name, critical_pressure, critical_temperature, accentric_factor, critical_compressibility, lennard_jonnes_diameter); }),
              py::arg("name"),
              py::arg("critical_pressure"),
              py::arg("critical_temperature"),
              py::arg("accentric_factor"),
              py::arg("critical_compressibility"),
              py::arg("lennard_jonnes_diameter"))
         .def(py::init([](std::string name,
                          double critical_pressure,
                          double critical_temperature,
                          double accentric_factor)
                       { return new Fluid(name, critical_pressure, critical_temperature, accentric_factor); }),
              py::arg("name"),
              py::arg("critical_pressure"),
              py::arg("critical_temperature"),
              py::arg("accentric_factor"))
         .def(py::init<>());

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
         .def("set_adsorbent", &PurePTA::SetAdsorbent, py::arg("Adsorbent"))
         .def("get_loading", &PurePTA::GetLoading,
              "Get single-component loading",
              py::arg("pressure"),
              py::arg("temperature"),
              py::arg("adsorption_potential_parameters"),
              py::arg("fluid"))
         .def("get_pressure", &PurePTA::GetPressure,
              "Get single-component pressure for a specific loading",
              py::arg("loading"),
              py::arg("temperature"),
              py::arg("adsorption_potential_parameters"),
              py::arg("fluid"),
              py::arg("P_estimate_"))
         .def("get_deviation_range", &PurePTA::GetDeviationRange,
              "Get absolute difference from experimental and calculated loadings",
              py::arg("deviation_type"),
              py::arg("experimental_loading"),
              py::arg("experimental_pressure"),
              py::arg("temperature"),
              py::arg("adsorption_potential_parameters"),
              py::arg("fluid"))
         .def("get_loadings", &PurePTA::GetLoadings,
              "Get single-component loading for multiple pressure datapoints",
              py::arg("experimental_pressure"),
              py::arg("temperature"),
              py::arg("adsorption_potential_parameters"),
              py::arg("fluid"));

     py::class_<MixturePTA>(m, "MixturePTA")
         .def(py::init<std::string &, std::string &, std::string &, std::size_t &>(),
              py::arg("adsorption_potential"),
              py::arg("equation_of_state"),
              py::arg("isotherm_type"),
              py::arg("number_of_layers"))
         .def_readwrite("adsorption_potential", &MixturePTA::potential)
         .def_readwrite("equation_of_state", &MixturePTA::equation_of_state)
         .def_readwrite("isotherm_type", &MixturePTA::isotherm_type)
         .def_readwrite("number_of_layers", &MixturePTA::num_of_layers)
         .def("get_loading", &MixturePTA::GetLoading,
              "Get multiple-component loadings",
              py::arg("composition"),
              py::arg("pressure"),
              py::arg("temperature"),
              py::arg("adsorption_potential_parameters"),
              py::arg("fluid_parameters"))
         .def("set_adsorbent", &MixturePTA::SetAdsorbent, py::arg("Adsorbent"));

     m.attr("DRA_POTENTIAL") = py::str(DRA_POTENTIAL);
     m.attr("STEELE_POTENTIAL") = py::str(STEELE_POTENTIAL);
     m.attr("LEE_POTENTIAL") = py::str(LEE_POTENTIAL);

#ifdef VERSION_INFO
     m.attr("__version__") = "0.0.1";
#else
     m.attr("__version__") = "dev";
#endif
}
