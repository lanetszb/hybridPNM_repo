/* MIT License
 *
 * Copyright (c) 2022 Aleksandr Zhuravlyov and Zakhar Lanets
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef HYBRIDPNM_WRAPPER_H
#define HYBRIDPNM_WRAPPER_H


#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "Props.h"
#include "Variables.h"
#include "loc/Loc.h"
#include "loc/LocP.h"
#include "loc/LocS.h"
#include "con/Con.h"
#include "con/ConP.h"
#include "con/ConS.h"
#include "equation/Equation.h"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace hybridPNM;

void modKFactor(Eigen::Ref<Eigen::VectorXd> KFactor, Eigen::Ref<Eigen::VectorXd> P) {
  KFactor.setIdentity();
  KFactor += P / 1000.;
}

PYBIND11_MODULE(hybridPNM, m) {

  py::class_<Props, std::shared_ptr<Props>>(m, "Props")
      .def(py::init<const std::map<std::string, std::variant<int, double, std::string>> &>(), "dict"_a)
      .def(py::init<>())
      .def("__repr__", &Props::print)
      .def_readwrite("_dict", &Props::_dict);


  py::class_<Variables, std::shared_ptr<Variables>>(m, "Variables")
      .def(py::init<const int &, const int &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "iPrev"_a, "iCurr"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a)

      .def_readwrite("_iPrev", &Variables::_iPrev)
      .def_readwrite("_iCurr", &Variables::_iCurr)
      .def_property("u", &Variables::getU, &Variables::setU)
      .def_property("dS0", &Variables::getDS0, &Variables::setDS0)
      .def_property("P", &Variables::getP, &Variables::setP)
      .def_property("S0", &Variables::getS0, &Variables::setS0)
      .def_property("data", &Variables::getData, &Variables::setData);


  py::class_<Loc, std::shared_ptr<Loc>, Variables>(m, "Loc")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a)

      .def("preComputations", &Loc::preComputations)
      .def("postComputations", &Loc::postComputations)
      .def("iterateIndices", &Loc::iterateIndices)
      .def("calculateCoefficients", &Loc::calculateCoefficients)
      .def("imposeBcDirichlet", &Loc::imposeBcDirichlet)

      .def_readwrite("_dV", &Loc::_dV)
      .def_readwrite("_free", &Loc::_free)
      .def_readwrite("_element", &Loc::_element)
      .def_readwrite("_bc", &Loc::_bc);


  py::class_<LocP, std::shared_ptr<LocP>, Loc>(m, "LocP")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a);


  py::class_<LocS, std::shared_ptr<LocS>, Loc>(m, "LocS")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a);


  py::class_<Con, std::shared_ptr<Con>, Variables>(m, "Con")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a)

      .def("preComputations", &Con::preComputations)
      .def("postComputations", &Con::postComputations)
      .def("iterateIndices", &Con::iterateIndices)
      .def("calculateCoefficients", &Con::calculateCoefficients)
      .def("imposeBcNeumann", &Con::imposeBcNeumann)
      .def("zeroCoefficients", &Con::zeroCoefficients, "names"_a)
      .def("calculateDS0", &Con::calculateDS0)

      .def_readwrite("_dV", &Con::_dV)
      .def_readwrite("_dA", &Con::_dA)
      .def_readwrite("_dL", &Con::_dL)
      .def_readwrite("_free", &Con::_free)
      .def_readwrite("_elements", &Con::_elements)
      .def_readwrite("_bc", &Con::_bc);


  py::class_<ConP, std::shared_ptr<ConP>, Con>(m, "ConP")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a)

      .def("calculateVisc", &ConP::calculateVisc)
      .def("calculateU", &ConP::calculateU)

      .def_readwrite("_K", &ConP::_K)
      .def_readwrite("_visc", &ConP::_visc)
      .def_readwrite("_PCMag", &ConP::_PCMag);


  py::class_<ConS, std::shared_ptr<ConS>, Con>(m, "ConS")
      .def(py::init<netgrid::Netgrid &, Props &,
               Eigen::Ref<Eigen::VectorXd>,
               Eigen::Ref<Eigen::VectorXd>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::vector<Eigen::Ref<Eigen::VectorXd>>,
               std::map<std::string, Eigen::Ref<Eigen::VectorXd>>>(),
           "netgrid"_a, "properties"_a, "u"_a, "dS0"_a, "P"_a, "S0"_a, "data"_a)

      .def("interpolate", &ConS::interpolate)

      .def_readwrite("_upwindCell", &ConS::_upwindCell);


  py::class_<Equation, std::shared_ptr<Equation>, Variables>(m, "Equation")
      .def(py::init<Loc &, Con &, std::vector<Eigen::Ref<Eigen::VectorXd>>>(),
           "loc"_a, "con"_a, "X"_a)

      .def("preComputations", &Equation::preComputations)
      .def("postComputations", &Equation::postComputations)
      .def("iterateIndices", &Equation::iterateIndices)
      .def("calculateCoefficients", &Equation::calculateCoefficients)
      .def("fillMatrix", &Equation::fillMatrix)
      .def("solve", &Equation::solve)
      .def("calculateData", &Equation::calculateData)

      .def_property("X", &Equation::getX, &Equation::setX)
      .def_readwrite("_matrix", &Equation::_matrix)
      .def_readwrite("_free", &Equation::_free)
      .def_readwrite("_bcTypes", &Equation::_bcTypes)
      .def_readwrite("_bc", &Equation::_bc)
      .def_readwrite("_tolerance", &Equation::_tolerance)
      .def_readwrite("_solver", &Equation::_solver)
      .def_readonly("_solverList", &Equation::_solverList);

  m.def("modKFactor", &modKFactor, "KFactor"_a, "P"_a);

}


#endif //HYBRIDPNM_WRAPPER_H