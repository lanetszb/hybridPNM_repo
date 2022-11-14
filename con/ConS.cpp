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


#include "ConS.h"

using namespace hybridPNM;


template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


ConS::ConS(netgrid::Netgrid &netgrid, Props &props,
           Eigen::Map<Eigen::VectorXd> u,
           Eigen::Map<Eigen::VectorXd> dS0,
           std::vector<Eigen::Map<Eigen::VectorXd>> P,
           std::vector<Eigen::Map<Eigen::VectorXd>> S0,
           std::map<std::string, Eigen::Map<Eigen::VectorXd>> data) :
    Con(netgrid, props, u, dS0, P, S0, data) {}

ConS::ConS(netgrid::Netgrid &netgrid, Props &props,
           Eigen::Ref<Eigen::VectorXd> u,
           Eigen::Ref<Eigen::VectorXd> dS0,
           std::vector<Eigen::Ref<Eigen::VectorXd>> P,
           std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
           std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data) :
    ConS(netgrid, props,
         setter<Eigen::VectorXd>(u),
         setter<Eigen::VectorXd>(dS0),
         vectorSetter<Eigen::VectorXd>(P),
         vectorSetter<Eigen::VectorXd>(S0),
         mapSetter<std::string, Eigen::VectorXd>(data)) {}

ConS::~ConS() {}


void ConS::preComputations() {
  interpolate();
}


void ConS::postComputations() {}


void ConS::calculateCoefficients() {

  auto xi = std::get<double>(_props._dict["xi"]);

  for (int cell = 0; cell < _netgrid._cellsN; cell++) {
    _free[cell] = 0;
    for (auto &[cellNeighbor, value]: _elements[cell])
      value = 0;
  }

  auto &S0 = _S0[_iPrev];

  for (auto &[face, cell]: _upwindCell) {

    auto itfFactor = xi * _u[face] * _dA[face];

    auto normal = _netgrid._normalsNeighborsCells[face][0];

    if (_netgrid._neighborsCells[face].size() == 1) {

      _elements[cell][cell] += normal * _u[face] * _dA[face];


      // implicit interface compression
      _elements[cell][cell] -= normal * itfFactor * (1. - 2. * S0[cell]);
      _free[cell] -= normal * itfFactor * S0[cell] * S0[cell];

      // explicit interface compression
      // _free[cell] -= normal * itfFactor * S0[cell] * (1. - S0[cell]);

    } else if (_netgrid._neighborsCells[face].size() == 2) {
      auto cellNeighbor = _netgrid._neighborsCells[face][1];
      auto normalNeighbor = _netgrid._normalsNeighborsCells[face][1];

      if (cell != _netgrid._neighborsCells[face][0]) {
        normal = _netgrid._normalsNeighborsCells[face][1];
        cellNeighbor = _netgrid._neighborsCells[face][0];
        normalNeighbor = _netgrid._normalsNeighborsCells[face][0];
      }

      _elements[cell][cell] += normal * _u[face] * _dA[face];
      _elements[cellNeighbor][cell] += normalNeighbor * _u[face] * _dA[face];

      // implicit interface compression
      _elements[cell][cell] -= normal * itfFactor * (1. - 2. * S0[cell]);
      _elements[cellNeighbor][cell] -= normalNeighbor * itfFactor * (1. - 2. * S0[cell]);
      _free[cell] -= normal * itfFactor * S0[cell] * S0[cell];
      _free[cellNeighbor] -= normalNeighbor * itfFactor * S0[cell] * S0[cell];

      // explicit interface compression
      // _free[cell] -= normal * itfFactor * S0[cell] * (1. - S0[cell]);
      // _free[cellNeighbor] -= normalNeighbor * itfFactor * S0[cell] * (1. - S0[cell]);

    }

  }

}


void ConS::imposeBcNeumann() {
  auto xi = std::get<double>(_props._dict["xi"]);
  auto &S0 = _S0[_iPrev];
  // value here acts as saturation!
  for (auto &[name, value]: _bc)
    for (auto &face: _netgrid._typesFaces[name])
      if (_u[face] > 0) {

        auto normal = _netgrid._normalsNeighborsCells[face][0];
        auto cell = _netgrid._neighborsCells[face][0];

        _free[cell] += normal * value * _u[face] * _dA[face];

        auto itfFactor = xi * _u[face] * _dA[face];

        // implicit interface compression
        _elements[cell][cell] -= normal * itfFactor * (1. - 2. * S0[cell]);
        _free[cell] -= normal * itfFactor * S0[cell] * S0[cell];

        // explicit interface compression
        // _free[cell] -= normal * itfFactor * S0[cell] * (1. - S0[cell]);

      }
}


void ConS::interpolate() {

  _upwindCell.clear();

  for (auto &face: _netgrid._typesFaces["nonbound"]) {
    auto &cell0 = _netgrid._neighborsCells[face][0];
    auto &cell1 = _netgrid._neighborsCells[face][1];
    auto &normal0 = _netgrid._normalsNeighborsCells[face][0];
    auto &normal1 = _netgrid._normalsNeighborsCells[face][1];

    if (normal0 * _u[face] > 0)
      _upwindCell[face] = cell0;
    else if (normal0 * _u[face] < 0)
      _upwindCell[face] = cell1;
  }

  for (auto &face: _netgrid._typesFaces["bound"]) {
    auto &cell = _netgrid._neighborsCells[face][0];
    auto &normal = _netgrid._normalsNeighborsCells[face][0];
    if (_u[face] < 0)
      _upwindCell[face] = cell;
  }

}