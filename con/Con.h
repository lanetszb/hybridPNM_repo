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
 @brief File for Con class.
 It is header which contains Con class.
*/

#ifndef HYBRIDPNM_CON_H
#define HYBRIDPNM_CON_H

#include <string>
#include <vector>

#include "Variables.h"
#include <Props.h>
#include <Netgrid.h>

namespace hybridPNM {

  /// This class is realisation of convection (surface) terms for elliptic partial differential equation.
  class Con : public Variables {

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
    Con(netgrid::Netgrid &netgrid, Props &props,
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
    Con(netgrid::Netgrid &netgrid, Props &props,
        Eigen::Ref<Eigen::VectorXd> u,
        Eigen::Ref<Eigen::VectorXd> dS0,
        std::vector<Eigen::Ref<Eigen::VectorXd>> P,
        std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
        std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data);

    /**
     Destructor is set by default.
    */
    virtual ~Con();


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
     Supposed to impose Neumann type boundary conditions.
    */
    virtual void imposeBcNeumann() {}

    /**
     Assign zeros to coefficients for matrix regarding particular face types.
     @param[in] names is types of faces.
    */
    void zeroCoefficients(const std::set<std::string> &names);

    /**
     Calculate delta of saturation of 0 phase regarding a face.
    */
    void calculateDS0();


    netgrid::Netgrid &_netgrid; ///< mesh.
    Props &_props; ///< Different properties.

    std::vector<double> _dV; ///< Cell volume.
    std::vector<double> _dA; ///< Face area.
    std::vector<double> _dL; ///< Cell centers distance with respect to face.

    std::vector<double> _free; ///< Free term coefficient.
    std::vector<std::map<uint32_t, double>> _elements; ///< Matrix coefficients.

    std::map<std::string, double> _bc; ///< Boundary conditions.

  };

}

#endif //HYBRIDPNM_CON_H
