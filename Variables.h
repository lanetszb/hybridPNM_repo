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
 @brief File for Variables class.
 It is header which contains Variables class.
*/


#ifndef HYBRIDPNM_VARIABLES_H
#define HYBRIDPNM_VARIABLES_H


#include <vector>
#include <Eigen/Dense>

#include "eigenSetGet.h"


namespace hybridPNM {

  /// This class is realisation of general attributes and it's accessors and mutators.
  class Variables {
  public:

    /**
     Constructor sets iPrev, iCurr, u, dS0, P, S0 and data
     where eigen vectors type is Map.
     @param[in] iPrev is previous time index.
     @param[in] iCurr is current time index.
     @param[in] u is velocity.
     @param[in] dS0 is delta of saturation of 0 phase regarding a face.
     @param[in] P is pressure.
     @param[in] S0 is saturation of 0 phase.
     @param[in] data is structure for interchange between modules.
    */
    Variables(const int &iPrev, const int &iCurr,
              Eigen::Map<Eigen::VectorXd> u,
              Eigen::Map<Eigen::VectorXd> dS0,
              std::vector<Eigen::Map<Eigen::VectorXd>> P,
              std::vector<Eigen::Map<Eigen::VectorXd>> S0,
              std::map<std::string, Eigen::Map<Eigen::VectorXd>> data);

    /**
     Constructor sets iPrev, iCurr, u, dS0, P, S0 and data
     where eigen vectors type is Ref.
     @param[in] iPrev is previous time index.
     @param[in] iCurr is current time index.
     @param[in] u is velocity.
     @param[in] dS0 is delta of saturation of 0 phase regarding a face.
     @param[in] P is pressure.
     @param[in] S0 is saturation of 0 phase.
     @param[in] data is structure for interchange between modules.
    */
    Variables(const int &iPrev, const int &iCurr,
              Eigen::Ref<Eigen::VectorXd> u,
              Eigen::Ref<Eigen::VectorXd> dS0,
              std::vector<Eigen::Ref<Eigen::VectorXd>> P,
              std::vector<Eigen::Ref<Eigen::VectorXd>> S0,
              std::map<std::string, Eigen::Ref<Eigen::VectorXd>> data);

    virtual ~Variables();

    int _iPrev; ///< Previous time index.
    int _iCurr; ///< Current time index.

    Eigen::Map<Eigen::VectorXd> _u; ///< Velocity.
    Eigen::Map<Eigen::VectorXd> _dS0; ///< Delta of saturation of 0 phase regarding a face.
    std::vector<Eigen::Map<Eigen::VectorXd>> _P; ///< Pressure.
    std::vector<Eigen::Map<Eigen::VectorXd>> _S0; ///< Saturation of 0 phase.

    std::map<std::string, Eigen::Map<Eigen::VectorXd>> _data; ///< Structure for interchange between modules.


    /**
     Accessor for velocity.
    */
    Eigen::Ref<Eigen::VectorXd> getU();

    /**
     Mutators for velocity.
     @param[in] u is velocity.
    */
    void setU(Eigen::Ref<Eigen::VectorXd> u);


    /**
     Accessor for delta of saturation of 0 phase regarding a face.
    */
    Eigen::Ref<Eigen::VectorXd> getDS0();

    /**
     Mutators for delta of saturation of 0 phase regarding a face.
     @param[in] dS0 is delta of saturation of 0 phase regarding a face.
    */
    void setDS0(Eigen::Ref<Eigen::VectorXd> dS0);


    /**
     Accessor for pressure.
    */
    std::vector<Eigen::Ref<Eigen::VectorXd>> getP();

    /**
     Mutators for pressure.
     @param[in] P is pressure.
    */
    void setP(std::vector<Eigen::Ref<Eigen::VectorXd>> &P);


    /**
     Accessor for saturation of 0 phase.
    */
    std::vector<Eigen::Ref<Eigen::VectorXd>> getS0();

    /**
     Mutators for saturation of 0 phase.
     @param[in] S0 is saturation of 0 phase.
    */
    void setS0(std::vector<Eigen::Ref<Eigen::VectorXd>> &S0);


    /**
     Accessor for structure for interchange between modules.
    */
    std::map<std::string, Eigen::Ref<Eigen::VectorXd>> getData();

    /**
     Mutators for structure for interchange between modules.
     @param[in] data structure for interchange between modules.
    */
    void setData(std::map<std::string, Eigen::Ref<Eigen::VectorXd>> &data);

  };

}

#endif //HYBRIDPNM_VARIABLES_H
