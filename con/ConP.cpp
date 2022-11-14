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


#include "ConP.h"


using namespace hybridPNM;


template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


ConP::ConP(netgrid::Netgrid &netgrid, Props &props,
           Eigen::Map<Eigen::VectorXd> u,
           Eigen::Map<Eigen::VectorXd> dS0,
           std::vector<Eigen::Map<Eigen::VectorXd>> P,
           std::vector<Eigen::Map<Eigen::VectorXd>> S0,
           std::map<std::string, Eigen::Map<Eigen::VectorXd>> data) :
    Con(netgrid, props, u, dS0, P, S0, data),
    _visc(std::vector<double>(_netgrid._cellsN, 0)),
    _K(std::vector<double>(_netgrid._cellsN, 0)),
    _PCMag(std::vector<double>(_netgrid._facesN, 0)) {

  auto dim = std::get<std::string>(_props._dict["dim"]);
  if (dim != "2D" and dim != "3D")
    throw std::invalid_argument("Received incorrect dim. Should be 2D or 3D.");

  /// ToDo: make sure that 3D _K is calculated correct
  for (auto &[throat, cells]: _netgrid._throatsCells) {
    auto &width = _netgrid._throatsWidths[throat];
    auto &depth = _netgrid._throatsDepths[throat];
    for (auto &cell: cells)
      if (dim == "2D")
        _K[cell] = width * width / 12.;
      else if (dim == "3D")
        _K[cell] = 5. * width * width * depth * depth / 18. / (width * width + depth * depth);
  }

  auto theta = std::get<double>(_props._dict["theta"]);
  auto beta = std::get<double>(_props._dict["beta"]);
  for (auto &[throat, faces]: _netgrid._throatsFaces) {
    auto &width = _netgrid._throatsWidths[throat];
    auto &depth = _netgrid._throatsDepths[throat];
    for (auto &face: faces) {
      auto factor = 2. * beta * cos(theta);
      if (dim == "2D")
        _PCMag[face] = factor * (1. / width);
      else if (dim == "3D")
        _PCMag[face] = factor * (1. / width + 1. / depth);
    }
  }

}


ConP::ConP(netgrid::Netgrid &netgrid, Props &props,
           Eigen::Ref<Eigen::VectorXd> u,
           Eigen::Ref<Eigen::VectorXd> dS0,
           std::vector<Eigen::Ref<Eigen::VectorXd>> P,
           std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
           std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data) :
    ConP(netgrid, props,
         setter<Eigen::VectorXd>(u),
         setter<Eigen::VectorXd>(dS0),
         vectorSetter<Eigen::VectorXd>(P),
         vectorSetter<Eigen::VectorXd>(S0),
         mapSetter<std::string, Eigen::VectorXd>(data)) {}


ConP::~ConP() {}


void ConP::preComputations() {
  calculateVisc();
  calculateDS0();
}

void ConP::postComputations() {
  calculateU();
}


void ConP::calculateCoefficients() {

  for (int cell = 0; cell < _netgrid._cellsN; cell++) {
    _free[cell] = 0;
    for (auto &[cellNeighbor, value]: _elements[cell])
      value = 0;
  }

  auto gamma = std::get<double>(_props._dict["gamma"]);

  for (auto &face: _netgrid._typesFaces["nonbound"]) {

    auto &cell0 = _netgrid._neighborsCells[face][0];
    auto &cell1 = _netgrid._neighborsCells[face][1];
    auto &normal0 = _netgrid._normalsNeighborsCells[face][0];
    auto &normal1 = _netgrid._normalsNeighborsCells[face][1];

    auto factor = conductanceAv(_data.at("KFactor")[cell0] * _K[cell0] / _visc[cell0],
                                _data.at("KFactor")[cell1] * _K[cell1] / _visc[cell1]) * _dA[face] / _dL[face];

    auto PC = factor * sgn(_dS0[face]) * std::pow(fabs(_dS0[face]), gamma) * _PCMag[face];

    _free[cell0] += normal0 * PC;
    _free[cell1] += normal1 * PC;

    _elements[cell0][cell0] += normal0 * normal0 * factor;
    _elements[cell0][cell1] += normal0 * normal1 * factor;

    _elements[cell1][cell1] += normal1 * normal1 * factor;
    _elements[cell1][cell0] += normal1 * normal0 * factor;

  }

}


void ConP::imposeBcNeumann() {

  for (auto &[name, value]: _bc)
    for (auto &face: _netgrid._typesFaces[name]) {
      auto cell = _netgrid._neighborsCells[face][0];
      auto normal = _netgrid._normalsNeighborsCells[face][0];
      _free[cell] += normal * value;
    }
}


void ConP::calculateVisc() {
  auto visc0 = std::get<double>(_props._dict["visc0"]);
  auto visc1 = std::get<double>(_props._dict["visc1"]);
  for (int cell = 0; cell < _netgrid._cellsN; cell++)
    _visc[cell] = _S0[_iPrev][cell] * visc0 + (1. - _S0[_iPrev][cell]) * visc1;
}


double ConP::conductanceAv(const double &value0, const double &value1) {
  return 2. * value0 * value1 / (value0 + value1);
}


void ConP::calculateU() {

  auto gamma = std::get<double>(_props._dict["gamma"]);

  for (auto &face: _netgrid._typesFaces["nonbound"]) {
    auto &cell0 = _netgrid._neighborsCells[face][0];
    auto &cell1 = _netgrid._neighborsCells[face][1];
    auto &normal0 = _netgrid._normalsNeighborsCells[face][0];
    auto &normal1 = _netgrid._normalsNeighborsCells[face][1];

    auto factor = conductanceAv(_data.at("KFactor")[cell0] * _K[cell0] / _visc[cell0],
                                _data.at("KFactor")[cell1] * _K[cell1] / _visc[cell1]) / _dL[face];

    auto PC = factor * sgn(_dS0[face]) * std::pow(fabs(_dS0[face]), gamma) * _PCMag[face];

    _u[face] = factor * (normal0 * _P[_iCurr][cell0] + normal1 * _P[_iCurr][cell1]) + PC;

  }

  for (auto &face: _netgrid._typesFaces["bound"]) {

    auto &cellNeighbor = _netgrid._neighborsCells[face][0];
    std::set < uint32_t > faceNeighborIdxs;
    uint32_t faceIdx;
    auto &faces = _netgrid._neighborsFaces[cellNeighbor];
    auto &normals = _netgrid._normalsNeighborsFaces[cellNeighbor];
    for (int idx = 0; idx < faces.size(); idx++)
      if (faces[idx] != face)
        faceNeighborIdxs.insert(idx);
      else
        faceIdx = idx;

    _u[face] = 0;
    for (auto &idx: faceNeighborIdxs)
      _u[face] -= _u[faces[idx]] * _dA[faces[idx]] / _dA[faces[faceIdx]] *
                  normals[idx] / normals[faceIdx];

  }

  for (auto &face: _netgrid._typesFaces["deadend"])
    _u[face] = 0;

  for (auto &[name, value]: _bc)
    for (auto &face: _netgrid._typesFaces[name]) {
      auto normal = _netgrid._normalsNeighborsCells[face][0];
      _u[face] = -normal * value / _dA[face];
    }

  std::set < std::string > namesBc{"inlet", "outlet"};
  for (auto &name: namesBc)
    if (not _bc.contains(name))
      for (auto &face: _netgrid._typesFaces[name]) {
        auto &cell = _netgrid._neighborsCells[face][0];
        auto &normal = _netgrid._normalsNeighborsCells[face][0];
        /// ToDo: make sure that sign below is correct
        _u[face] += normal * _data.at("qM")[cell] / _dA[face];
      }

}