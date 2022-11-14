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
 @brief File for Props class.
 It is header which contains Props class.
*/

#ifndef HYBRIDPNM_PROPS_H
#define HYBRIDPNM_PROPS_H

#include <iostream>
#include <map>
#include <variant>
#include <string>

namespace hybridPNM {

  /// This class is realisation of different properties.
  class Props {

  public:

    /**
     Constructor sets by default.
    */
    Props();

    /**
     Constructor sets dictionary with properties.
     @param[in] dict contains different properties.
    */
    Props(const std::map<std::string, std::variant<int, double, std::string>> &dict);

    /**
     Destructor is set by default.
    */
    virtual ~Props();

    /**
     Print properties.
     @return string with properties.
    */
    std::string print() const;

    std::map<std::string, std::variant<int, double, std::string>> _dict; ///< Dictionary with properties.

  };

}

/**
 Overload insertion operator.
 @param[in] stream is the instance of class std::ostream.
 @param[in] properties is the instance of this class.
 @return instance of class std::ostream..
*/
std::ostream &operator<<(std::ostream &stream, const hybridPNM::Props &properties);


#endif //HYBRIDPNM_PROPS_H


