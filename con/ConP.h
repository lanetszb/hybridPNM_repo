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
 @brief File for ConP class.
 It is header which contains ConP class.
*/

#ifndef HYBRIDPNM_CONP_H
#define HYBRIDPNM_CONP_H

#include "Con.h"

namespace hybridPNM {

  /// This class is realisation of convection (surface) terms for pressure elliptic partial differential equation.
  class ConP : public Con {

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
    ConP(netgrid::Netgrid &netgrid, Props &props,
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
    ConP(netgrid::Netgrid &netgrid, Props &props,
         Eigen::Ref<Eigen::VectorXd> u,
         Eigen::Ref<Eigen::VectorXd> dS0,
         std::vector<Eigen::Ref<Eigen::VectorXd>> P,
         std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
         std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data);

    /**
     Destructor is set by default.
    */
    virtual ~ConP();


    /**
     Make for different precomputations.
    */
    virtual void preComputations() override final;

    /**
     Make for different postcomputations.
    */
    virtual void postComputations() override final;

    /**
     Calculate coefficients for target SLA.
    */
    virtual void calculateCoefficients() override final;

    /**
     Impose Neumann type boundary conditions.
    */
    virtual void imposeBcNeumann() override final;

    /**
     Calculate viscosity.
    */
    void calculateVisc();

    /**
     Calculate average conductance.
     @param[in] value0 conductance value 0.
     @param[in] value1 conductance value 1.
     @return average conductance.
    */
    double conductanceAv(const double &value0, const double &value1);

    /**
     Calculate velocity.
    */
    void calculateU();


    std::vector<double> _visc; ///< Viscosity.

    std::vector<double> _K; ///< Permeability.


    std::vector<double> _PCMag; ///< Capillary pressure magnitude.

  };

}

#endif //HYBRIDPNM_CONP_H
