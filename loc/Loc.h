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
 @brief File for Loc class.
 It is header which contains Loc class.
*/


#ifndef HYBRIDPNM_LOC_H
#define HYBRIDPNM_LOC_H


#include "Variables.h"
#include "Props.h"
#include <Netgrid.h>


namespace hybridPNM {

  /// This class is realisation of location (volume) terms for elliptic partial differential equation.
  class Loc : public Variables {

  public:

    /**
     Constructor sets netgrid, props, u, dS0, P, S0 and data
     where eigen vectors type is Map.
     @param[in] netgrid is mesh.
     @param[in] props is different properties.
     @param[in] u is velocity.
     @param[in] dS0 is delta of saturation of 0 phase regarding a face.
     @param[in] P is pressure.
     @param[in] S0 is saturation of 0 phase.
     @param[in] data is structure for interchange between modules.
    */
    Loc(netgrid::Netgrid &netgrid, Props &props,
        Eigen::Map<Eigen::VectorXd> u,
        Eigen::Map<Eigen::VectorXd> dS0,
        std::vector<Eigen::Map<Eigen::VectorXd>> P,
        std::vector<Eigen::Map<Eigen::VectorXd>> S0,
        std::map<std::string, Eigen::Map<Eigen::VectorXd>> data);

    /**
     Constructor sets netgrid, props, u, dS0, P, S0 and data
     where eigen vectors type is Ref.
     @param[in] netgrid is mesh.
     @param[in] props is different properties.
     @param[in] u is velocity.
     @param[in] dS0 is delta of saturation of 0 phase regarding a face.
     @param[in] P is pressure.
     @param[in] S0 is saturation of 0 phase.
     @param[in] data is structure for interchange between modules.
    */
    Loc(netgrid::Netgrid &netgrid, Props &props,
        Eigen::Ref<Eigen::VectorXd> u,
        Eigen::Ref<Eigen::VectorXd> dS0,
        std::vector<Eigen::Ref<Eigen::VectorXd>> P,
        std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
        std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data);

    /**
     Destructor is set by default.
    */
    virtual ~Loc();

    /**
     Supposed to be responsible for different precomputations.
    */
    virtual void preComputations() {}

    /**
     Supposed to be responsible for different postcomputations.
    */
    virtual void postComputations() {}

    /**
     Iterate time indices (0,1).
    */
    void iterateIndices();

    /**
     Supposed to calculate coefficients for target SLA.
    */
    virtual void calculateCoefficients() {}

    /**
     Impose Dirichlet type boundary conditions.
    */
    void imposeBcDirichlet();

    netgrid::Netgrid &_netgrid; ///< mesh.
    Props &_props; ///< Different properties.

    std::vector<double> _dV; ///< Cell volume.

    std::vector<double> _free; ///< Free term coefficient.
    std::vector<double> _element; ///< Main diagonal coefficient.

    std::map<std::string, double> _bc; ///< Boundary conditions.

  };

}

#endif //HYBRIDPNM_LOC_H
