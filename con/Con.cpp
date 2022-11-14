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


#include "Con.h"
#include "eigenSetGet.h"


using namespace hybridPNM;


Con::Con(netgrid::Netgrid &netgrid, Props &props,
         Eigen::Map<Eigen::VectorXd> u,
         Eigen::Map<Eigen::VectorXd> dS0,
         std::vector<Eigen::Map<Eigen::VectorXd>> P,
         std::vector<Eigen::Map<Eigen::VectorXd>> S0,
         std::map<std::string, Eigen::Map<Eigen::VectorXd>> data)
    : Variables(0, 1, u, dS0, P, S0, data), _netgrid(netgrid), _props(props),
      _dV(std::vector<double>(_netgrid._cellsN, 0)),
      _dA(std::vector<double>(_netgrid._facesN, 0)),
      _dL(std::vector<double>(_netgrid._facesN, 0)),
      _free(std::vector<double>(_netgrid._cellsN, 0)),
      _elements(_netgrid._cellsN) {

  for (auto &[throat, cells]: _netgrid._throatsCells)
    for (auto &cell: cells)
      _dV[cell] = _netgrid._throatsDVs[throat];

  for (auto &[throat, faces]: _netgrid._throatsFaces)
    for (auto &face: faces) {
      _dA[face] = _netgrid._throatsAs[throat];
      _dL[face] = _netgrid._throatsDLs[throat];
    }

  for (auto &[cell, faces]: _netgrid._neighborsFaces)
    for (auto &face: faces)
      for (auto &cellNeighbor: _netgrid._neighborsCells[face])
        _elements[cell][cellNeighbor] = 0;

}


Con::Con(netgrid::Netgrid &netgrid, Props &properties,
         Eigen::Ref<Eigen::VectorXd> u,
         Eigen::Ref<Eigen::VectorXd> dS0,
         std::vector<Eigen::Ref<Eigen::VectorXd>> P,
         std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
         std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data) :
    Con(netgrid, properties,
        setter<Eigen::VectorXd>(u),
        setter<Eigen::VectorXd>(dS0),
        vectorSetter<Eigen::VectorXd>(P),
        vectorSetter<Eigen::VectorXd>(S0),
        mapSetter<std::string, Eigen::VectorXd>(data)) {}


Con::~Con() {}


void Con::iterateIndices() {
  std::swap(_iPrev, _iCurr);
}


void Con::zeroCoefficients(const std::set<std::string> &names) {

  for (auto &name: names)
    for (auto &face: _netgrid._typesFaces[name])
      for (auto &cell: _netgrid._neighborsCells[face]) {
        _free[cell] = 0;
        for (auto &[cellNeighbor, value]: _elements[cell])
          value = 0;
      }

}


void Con::calculateDS0() {

  for (auto &face: _netgrid._typesFaces["nonbound"]) {
    auto &cell0 = _netgrid._neighborsCells[face][0];
    auto &cell1 = _netgrid._neighborsCells[face][1];
    auto &normal0 = _netgrid._normalsNeighborsCells[face][0];
    auto &normal1 = _netgrid._normalsNeighborsCells[face][1];

    _dS0[face] = normal0 * _S0[_iPrev][cell0] + normal1 * _S0[_iPrev][cell1];

  }

  for (auto &face: _netgrid._typesFaces["bound"]) {
    auto &cellNeighbor = _netgrid._neighborsCells[face][0];
    std::set<uint32_t> faceNeighborIdxs;
    uint32_t faceIdx;
    auto &faces = _netgrid._neighborsFaces[cellNeighbor];
    auto &normals = _netgrid._normalsNeighborsFaces[cellNeighbor];
    for (int idx = 0; idx < faces.size(); idx++)
      if (faces[idx] != face)
        faceNeighborIdxs.insert(idx);
      else
        faceIdx = idx;

    double QTotal = 0;
    for (auto &idx: faceNeighborIdxs)
      QTotal += fabs(_u[faces[idx]] * _dA[faces[idx]]);

    _dS0[face] = 0;
    if (QTotal > 0)
      for (auto &idx: faceNeighborIdxs)
        _dS0[face] -= _dS0[faces[idx]] * normals[idx] / normals[faceIdx] *
                      fabs(_u[faces[idx]]) * _dA[faces[idx]] / QTotal;
    else
      for (auto &idx: faceNeighborIdxs)
        _dS0[face] -= _dS0[faces[idx]] * normals[idx] / normals[faceIdx]
                      / double(faceNeighborIdxs.size());
  }

}