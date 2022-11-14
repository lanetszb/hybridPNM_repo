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


#include "Props.h"

#include <sstream>


using namespace hybridPNM;


Props::Props() {};

Props::Props(const std::map<std::string, std::variant<int, double, std::string>> &dict) : _dict(dict) {}

Props::~Props() {}


std::string Props::print() const {

  std::stringstream stream;
  stream << "{";
  for (auto &ent: _dict) {
    stream << "'" << ent.first << "'" << ": ";
    if (std::get_if<int>(&ent.second))
      stream << std::get<int>(ent.second) << ", ";
    else if (std::get_if<double>(&ent.second))
      stream << std::get<double>(ent.second) << ", ";
    else if (std::get_if<std::string>(&ent.second))
      stream << "'" << std::get<std::string>(ent.second) << "'" << ", ";
  }
  stream.seekp(-2, std::ios_base::end);
  stream << "}";

  return stream.str();
}


std::ostream &operator<<(std::ostream &stream, const Props &properties) {
  stream << properties.print();
  return stream;
}


