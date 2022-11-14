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


#include "Variables.h"


using namespace hybridPNM;


Variables::Variables(const int &iPrev, const int &iCurr,
                     Eigen::Map<Eigen::VectorXd> u,
                     Eigen::Map<Eigen::VectorXd> dS0,
                     std::vector<Eigen::Map<Eigen::VectorXd>> P,
                     std::vector<Eigen::Map<Eigen::VectorXd>> S0,
                     std::map<std::string, Eigen::Map<Eigen::VectorXd>> data)
    : _iPrev(iPrev), _iCurr(iCurr), _u(u), _dS0(dS0), _P(P), _S0(S0), _data(data) {}


Variables::Variables(const int &iPrev, const int &iCurr,
                     Eigen::Ref<Eigen::VectorXd> u,
                     Eigen::Ref<Eigen::VectorXd> dS0,
                     std::vector<Eigen::Ref<Eigen::VectorXd>> P,
                     std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
                     std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data) :
    Variables(iPrev, iCurr,
              setter<Eigen::VectorXd>(u),
              setter<Eigen::VectorXd>(dS0),
              vectorSetter<Eigen::VectorXd>(P),
              vectorSetter<Eigen::VectorXd>(S0),
              mapSetter<std::string, Eigen::VectorXd>(data)) {}


Variables::~Variables() {}


Eigen::Ref<Eigen::VectorXd> Variables::getU() {
  return _u;
}


void Variables::setU(Eigen::Ref<Eigen::VectorXd> u) {
  setter<Eigen::VectorXd>(u, _u);
}


Eigen::Ref<Eigen::VectorXd> Variables::getDS0() {
  return _dS0;
}


void Variables::setDS0(Eigen::Ref<Eigen::VectorXd> dS0) {
  setter<Eigen::VectorXd>(dS0, _dS0);
}


std::vector<Eigen::Ref<Eigen::VectorXd>> Variables::getP() {
  return vectorGetter<Eigen::VectorXd>(_P);
}


void Variables::setP(std::vector<Eigen::Ref<Eigen::VectorXd>> &P0) {
  vectorSetter<Eigen::VectorXd>(P0, _P);
}


std::vector<Eigen::Ref<Eigen::VectorXd>> Variables::getS0() {
  return vectorGetter<Eigen::VectorXd>(_S0);
}


void Variables::setS0(std::vector<Eigen::Ref<Eigen::VectorXd>> &S0) {
  vectorSetter<Eigen::VectorXd>(S0, _S0);
}


std::map<std::string, Eigen::Ref<Eigen::VectorXd>> Variables::getData() {
  return mapGetter<std::string, Eigen::VectorXd>(_data);
}

void Variables::setData(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> &data) {
  mapSetter<std::string, Eigen::VectorXd>(data, _data);
}