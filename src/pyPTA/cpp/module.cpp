#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include "src/data_classes.h"
#include "src/pta_pure.h"
#include "src/pta_mixture.h"

namespace py = pybind11;

PYBIND11_MODULE(cPTA, m)
{
	m.doc() = "C/C++ methods to run the PTA model";

	py::class_<Adsorbent>(m, "Adsorbent")
		.def(py::init<double &, double &>())
		.def_readwrite("lj_diameter", &Adsorbent::sigma_ss)
		.def_readwrite("atomic_density", &Adsorbent::rho_atoms);

	py::class_<Fluid>(m, "Fluid")
		.def(py::init<>())
		.def_readwrite("critical_pressure", &Fluid::critical_pressure)
		.def_readwrite("critical_temperature", &Fluid::critical_temperature)
		.def_readwrite("critical_compressibility", &Fluid::critical_compressibility)
		.def_readwrite("acentric_factor", &Fluid::acentric_factor)
		.def_readwrite("lj_diameter", &Fluid::lj_diameter)
		.def_readwrite("lj_energy", &Fluid::lj_energy);

	py::class_<PurePTA>(m, "PurePTA")
		.def(py::init<std::string &, std::string &, std::string &, std::size_t &>(),
			 py::arg("adsorption_potential"),
			 py::arg("equation_of_state"),
			 py::arg("isotherm_type"),
			 py::arg("number_of_layers"))
		.def_readwrite("adsorption_potential", &PurePTA::potential)
		.def_readwrite("equation_of_state", &PurePTA::equation_of_state)
		.def_readwrite("isotherm_type", &PurePTA::isotherm_type)
		.def_readwrite("number_of_layers", &PurePTA::num_of_layers)
		.def("get_loading", &PurePTA::get_loading,
			 "Get single-component loading",
			 py::arg("pressure"),
			 py::arg("temperature"),
			 py::arg("adsorption_potential_parameters"),
			 py::arg("fluid_parameters"))
		.def("get_absolute_deviation", &PurePTA::get_absolute_deviation,
			 "Get absolute difference from experimental and calculated loadings",
			 py::arg("deviation_type"),
			 py::arg("experimental_loading"),
			 py::arg("experimental_pressure"),
			 py::arg("temperature"),
			 py::arg("adsorption_potential_parameters"),
			 py::arg("fluid_parameters"))
		.def("get_multiple_loadings", &PurePTA::get_multiple_loadings,
			 "Get single-component loading for multiple pressure datapoints",
			 py::arg("experimental_pressure"),
			 py::arg("temperature"),
			 py::arg("adsorption_potential_parameters"),
			 py::arg("fluid_parameters"));

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
		.def("get_loading", &MixturePTA::get_loading,
			 "Get multiple-component loadings",
			 py::arg("composition"),
			 py::arg("pressure"),
			 py::arg("temperature"),
			 py::arg("adsorption_potential_parameters"),
			 py::arg("fluid_parameters"));
}