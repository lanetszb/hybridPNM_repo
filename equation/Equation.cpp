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


#include "Equation.h"

#include <set>


using namespace hybridPNM;


Equation::Equation(Loc &loc, Con &con, std::vector<Eigen::Map<Eigen::VectorXd>> X) :
    Variables(0, 1, loc._u, loc._dS0, loc._P, loc._S0, loc._data),
    _loc(loc),
    _con(con),
    _X(X),
    _matrix(_loc._netgrid._cellsN, _loc._netgrid._cellsN),
    _free(_loc._netgrid._cellsN),
    _tolerance(0),
    _solver("biCGSTAB"),
    _solverList({"biCGSTAB", "sparseLU"}) {

  auto cellsN = _loc._netgrid._cellsN;
  auto deadFacesN = _loc._netgrid._typesFaces["bound"].size();
  int nonneighborsConnectionsN = 0;
  for (auto &[cell, faces]: _loc._netgrid._neighborsFaces)
    if (faces.size() > 2)
      nonneighborsConnectionsN += faces.size() - 2;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * cellsN - deadFacesN + nonneighborsConnectionsN);
  for (auto &[cell, faces]: _loc._netgrid._neighborsFaces) {
    triplets.emplace_back(cell, cell);
    for (auto &face: faces)
      for (auto &cellNeighbor: _loc._netgrid._neighborsCells[face])
        if (cellNeighbor != cell)
          triplets.emplace_back(cell, cellNeighbor);
  }
  _matrix.setFromTriplets(triplets.begin(), triplets.end());

}


Equation::Equation(Loc &loc, Con &con, std::vector<Eigen::Ref<Eigen::VectorXd>> X) :
    Equation(loc, con, vectorSetter<Eigen::VectorXd>(X)) {}


Equation::~Equation() {}


void Equation::preComputations() {
  _loc.preComputations();
  _con.preComputations();
}

void Equation::postComputations() {
  _loc.postComputations();
  _con.postComputations();
}


void Equation::iterateIndices() {
  std::swap(_iPrev, _iCurr);
  _loc.iterateIndices();
  _con.iterateIndices();
}


void Equation::calculateCoefficients() {
  preComputations();
  _loc.calculateCoefficients();
  _con.calculateCoefficients();
  imposeBc();
}


void Equation::imposeBc() {

  std::map<std::string, double> _bcDirichlet;
  std::set<std::string> boundsDirichlet;
  std::map<std::string, double> _bcNeumann;
  for (const auto &[name, type]: _bcTypes)
    if (type == "Dirichlet") {
      _bcDirichlet[name] = _bc[name];
      boundsDirichlet.insert(name);
    } else if (type == "Neumann")
      _bcNeumann[name] = _bc[name];

  _loc._bc = _bcDirichlet;
  _loc.imposeBcDirichlet();
  _con.zeroCoefficients(boundsDirichlet);

  _con._bc = _bcNeumann;
  _con.imposeBcNeumann();

}


void Equation::fillMatrix() {

  for (int cell = 0; cell < _loc._netgrid._cellsN; cell++)
    _free[cell] = _loc._free[cell] + _con._free[cell];

  _matrix.setZero();

  for (int cell = 0; cell < _loc._netgrid._cellsN; cell++)
    _matrix.coeffRef(cell, cell) += _loc._element[cell];

  for (int cell = 0; cell < _loc._netgrid._cellsN; cell++)
    for (auto &[cellNeighbor, value]: _con._elements[cell])
      _matrix.coeffRef(cell, cellNeighbor) += value;

}


void Equation::solve() {

  if (_solver == "biCGSTAB") {

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> biCGSTAB;
    if (_tolerance != 0)
      biCGSTAB.setTolerance(_tolerance);
    biCGSTAB.compute(_matrix);
    _X[_iCurr] = biCGSTAB.solveWithGuess(-_free, _X[_iPrev]);

  } else if (_solver == "sparseLU") {

    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparseLU;
    sparseLU.compute(_matrix);
    _X[_iCurr] = sparseLU.solve(-_free);

  }

  postComputations();

}


void Equation::calculateData(){
  for (int cell = 0; cell < _loc._netgrid._cellsN; cell++){
    _data.at("P")[cell] = _P[_iCurr][cell];
    _data.at("S0")[cell] = _S0[_iCurr][cell];
  }
}


std::vector<Eigen::Ref<Eigen::VectorXd>> Equation::getX() {
  return vectorGetter<Eigen::VectorXd>(_X);
}


void Equation::setX(std::vector<Eigen::Ref<Eigen::VectorXd>> &X) {
  vectorSetter<Eigen::VectorXd>(X, _X);
}
