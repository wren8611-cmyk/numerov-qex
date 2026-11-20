// radialsolver.cxx
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>

#include "channel_system.hxx"
#include "potential_matrix.hxx"
#include "radial_solver.hxx"

#include <utility>
#include <sstream>

namespace py = pybind11;

PYBIND11_MODULE(radialsolver, m) {
    m.doc() = "Coupled-channel radial Schrodinger solver";

    py::enum_<BoundaryCondition>(m, "BoundaryCondition")
        .value("Outgoing", BoundaryCondition::Outgoing)
        .value("Incoming", BoundaryCondition::Incoming)
        .export_values();

    py::class_<ChannelInfo>(m, "ChannelInfo")
        .def(py::init<int, std::string, double, double, int, double, double>(),
             py::arg("id"), py::arg("name"), py::arg("m1"), py::arg("m2"),
             py::arg("ell") = 0, py::arg("z1") = 0.0, py::arg("z2") = 0.0)
        .def_readonly("id", &ChannelInfo::id)
        .def_readonly("name", &ChannelInfo::name)
        .def_readonly("mu", &ChannelInfo::mu)
        .def_readonly("threshold", &ChannelInfo::threshold)
        .def_readonly("ell", &ChannelInfo::ell)
        .def_readonly("charge_product", &ChannelInfo::charge_product);

    py::class_<ChannelSystem>(m, "ChannelSystem")
        .def(py::init<std::vector<ChannelInfo>>(), py::arg("channels"))
        .def("nchannels", &ChannelSystem::nchannels)
        .def("channel_index", &ChannelSystem::channel_index)
        .def("channel_by_id", &ChannelSystem::channel_by_id,
             py::return_value_policy::reference_internal);

    py::class_<PotentialInfo>(m, "PotentialInfo")
        .def(py::init([](std::string name, std::vector<double> pars,
                         int left, int right) {
            return PotentialInfo{std::move(name), std::move(pars), left, right};
        }), py::arg("name"), py::arg("pars"),
            py::arg("channelid_left"), py::arg("channelid_right"))
        .def_readonly("name", &PotentialInfo::name)
        .def_readonly("pars", &PotentialInfo::pars)
        .def_readonly("channelid_left", &PotentialInfo::channelid_left)
        .def_readonly("channelid_right", &PotentialInfo::channelid_right);

    py::class_<PotentialMatrix>(m, "PotentialMatrix")
        .def(py::init<const std::vector<PotentialInfo>&, const ChannelSystem&>(),
             py::arg("potentials"), py::arg("system"), py::keep_alive<1, 3>())
        .def("append", py::overload_cast<const PotentialInfo&>(&PotentialMatrix::append), py::arg("potential"))
        .def("append", py::overload_cast<const std::vector<PotentialInfo>&>(&PotentialMatrix::append), py::arg("potentials"))
        .def("clear", &PotentialMatrix::clear)
        .def("size", &PotentialMatrix::size)
        .def("display_info", [](const PotentialMatrix& self) {
            std::ostringstream oss;
            self.display_info(oss);
            py::print(oss.str(), py::arg("end") = ""); })
        .def("evaluate", &PotentialMatrix::evaluate, py::arg("r"));

    py::class_<MatrixFunction>(m, "MatrixFunction")
        .def("x", &MatrixFunction::x)
        .def("size", &MatrixFunction::size)
        .def("nrows", &MatrixFunction::nrows)
        .def("ncols", &MatrixFunction::ncols)
        .def("__len__", &MatrixFunction::size)
        .def("__getitem__", [](const MatrixFunction& wf, py::ssize_t index) {
            const auto n = static_cast<py::ssize_t>(wf.size());
            if (index < 0) index += n;
            if (index < 0 || index >= n) throw py::index_error();
            return Eigen::MatrixXcd(wf[static_cast<std::size_t>(index)]);
        });

    py::class_<ScatteringResult>(m, "ScatteringResult")
        .def_readonly("S", &ScatteringResult::S)
        .def_readonly("S_flux", &ScatteringResult::S_flux)
        .def_readonly("wavefunction", &ScatteringResult::wavefunction, py::return_value_policy::reference_internal)
        .def_property_readonly(
            "system",
            [](const ScatteringResult& self) -> const ChannelSystem& { return self.system; },
            py::return_value_policy::reference_internal
        )
        .def_readonly("E", &ScatteringResult::E)
        .def_readonly("open_channel_ids", &ScatteringResult::open_channel_ids)
        .def_readonly("boundary_condition", &ScatteringResult::boundary_condition)
        .def("open_channel_index", &ScatteringResult::open_channel_index, py::arg("channel_id"))
        .def("psi", &ScatteringResult::psi, py::arg("radial_index"), py::arg("channel_id"), py::arg("reference_channel_id"))
        .def("is_open", &ScatteringResult::is_open, py::arg("channel_id"))
        .def("s", &ScatteringResult::s, py::arg("channel_id_left"), py::arg("channel_id_right"))
        .def("s_flux", &ScatteringResult::s_flux, py::arg("channel_id_left"), py::arg("channel_id_right"));

    py::class_<RadialSolver>(m, "RadialSolver")
        .def(py::init<const ChannelSystem&, const PotentialMatrix&,
                      std::vector<double>, BoundaryCondition>(),
             py::arg("system"), py::arg("potential"), py::arg("rgrid"),
             py::arg("bc") = BoundaryCondition::Outgoing,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def("solve", &RadialSolver::solve, py::arg("E"),
             py::keep_alive<0, 1>());
}
