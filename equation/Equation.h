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

/**
 @file
 @brief File for Equation class.
 It is header which contains Equation class.
*/

#ifndef HYBRIDPNM_EQUATION_H
#define HYBRIDPNM_EQUATION_H


#include <Eigen/Sparse>

#include "Variables.h"
#include "con/Con.h"
#include "loc/Loc.h"

namespace hybridPNM {

  /// This class is realisation of elliptic partial differential equation.
  class Equation : public Variables {

  public:

    /**
     Constructor sets loc, con and X
     where eigen vectors type is Map.
     @param[in] loc is location (volume) term.
     @param[in] con is convection (surface) term.
     @param[in] X is equation variable.
    */
    Equation(Loc &loc, Con &con, std::vector<Eigen::Map<Eigen::VectorXd>> X);

    /**
     Constructor sets loc, con and X
     where eigen vectors type is Ref.
     @param[in] loc is location (volume) term.
     @param[in] con is convection (surface) term.
     @param[in] X is equation variable.
    */
    Equation(Loc &loc, Con &con, std::vector<Eigen::Ref<Eigen::VectorXd>> X);

    /**
     Destructor is set by default.
    */
    virtual ~Equation();

    /**
     Make for different precomputations.
    */
    void preComputations();

    /**
     Make for different postcomputations.
    */
    void postComputations();

    /**
     Iterate time indices (0,1).
    */
    void iterateIndices();

    /**
     Calculate coefficients for target SLA.
    */
    void calculateCoefficients();

    /**
     Impose boundary conditions.
    */
    void imposeBc();

    /**
     Fill target SLA matrix and free vector.
    */
    void fillMatrix();

    /**
     Solve target SLA.
    */
    void solve();

    /**
     Make calculations to process structure for interchange between modules.
    */
    void calculateData();


    Loc &_loc; ///< Location (volume) term.
    Con &_con; ///< Convection (surface) term.

    std::vector<Eigen::Map<Eigen::VectorXd>> _X; ///< Equation variable.

    Eigen::SparseMatrix<double, Eigen::RowMajor> _matrix; ///< Matrix for target SLE.
    Eigen::VectorXd _free; ///< Free term for target SLE.

    std::map<std::string, std::string> _bcTypes; ///< Boundary conditions types.
    std::map<std::string, double> _bc; ///< Boundary conditions.

    double _tolerance; ///< Tolerance for equation solution.

    std::string _solver; ///< Solver name.

    std::vector<std::string> _solverList; ///< List of available solvers.


    /**
     Accessor for equation variable.
    */
    std::vector<Eigen::Ref<Eigen::VectorXd>> getX();

    /**
     Mutators for equation variable.
     @param[in] X is equation variable.
    */
    void setX(std::vector<Eigen::Ref<Eigen::VectorXd>> &X);

  };

}

#endif //HYBRIDPNM_EQUATION_H
