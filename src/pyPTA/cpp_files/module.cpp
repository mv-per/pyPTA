#include "include/pybind11/pybind11.h"
#include "include/pybind11/embed.h"
#include "include/pybind11/stl.h"

#include "src/pta_pure.h"

namespace py = pybind11;

PYBIND11_MODULE(cpta, m)
{

	py::class_<PurePTA>(m, "PurePTA")
		.def(py::init<std::string &, std::string &, std::string &, std::size_t &>())
		.def("get_loading", &PurePTA::get_loading, "Get single-component loading")
		.def("get_absolute_deviation", &PurePTA::get_absolute_deviation, "Get absolute difference from experimental and calculated loadings")
		.def("get_multiple_loadings", &PurePTA::get_multiple_loadings, "Get single-component loading for multiple pressure datapoints");
}