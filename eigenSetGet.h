/* MIT License
 *
 * Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
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
 @brief File for accessors and mutators of Eigen attributes.
 It is header which contains accessors and mutators of Eigen attributes.
*/


#ifndef HYBRIDPNM_EIGENSETGET_H
#define HYBRIDPNM_EIGENSETGET_H

#include <Eigen/Dense>

#include <iostream>
#include <map>
#include <variant>


namespace hybridPNM {

  template<class T>
  void setter(Eigen::Ref<T> &source, Eigen::Map<T> &sink) {
    new(&sink) Eigen::Map<T>(source.data(), source.size());
  };

  template<class T>
  Eigen::Map<T> setter(Eigen::Ref<T> &source) {
    return Eigen::Map<T>(source.data(), source.size());
  };

  template<class T>
  std::vector<Eigen::Ref<T>>
  vectorGetter(std::vector<Eigen::Map<T>> &source) {

    std::vector<Eigen::Ref<T>> sink;
    for (auto value: source)
      sink.push_back(Eigen::Ref<T>(value));

    return sink;
  }

  template<class T>
  void vectorSetter(std::vector<Eigen::Ref<T>> &source,
                    std::vector<Eigen::Map<T>> &sink) {
    sink.clear();
    for (auto &value: source)
      sink.push_back(Eigen::Map<T>(value.data(), value.size()));
  }

  template<class T>
  std::vector<Eigen::Map<T>> vectorSetter(std::vector<Eigen::Ref<T>> &source) {
    std::vector<Eigen::Map<T>> sink;
    for (auto &value: source)
      sink.push_back(Eigen::Map<T>(value.data(), value.size()));
    return sink;
  }

  template<class K, class V>
  std::map<K, Eigen::Ref<V>> mapGetter(std::map<K, Eigen::Map<V>> &source) {
    std::map<K, Eigen::Ref<V>> sink;
    for (auto const &ent: source)
      sink.insert(std::pair<K, Eigen::Ref<V>>(ent.first, source.at(ent.first)));
    return sink;
  }

  template<class K, class V>
  void mapSetter(std::map<K, Eigen::Ref<V>> &source, std::map<K, Eigen::Map<V>> &sink) {
    sink.clear();
    for (auto const &ent: source)
      sink.insert(std::pair<K, Eigen::Map<V>>(ent.first, Eigen::Map<V>(source.at(ent.first).data(),
                                                                       source.at(ent.first).size())));
  }

  template<class K, class V>
  std::map<K, Eigen::Map<V>> mapSetter(std::map<K, Eigen::Ref<V>> &source) {
    std::map<K, Eigen::Map<V>> sink;
    for (auto const &ent: source)
      sink.insert(std::pair<K, Eigen::Map<V>>(ent.first, Eigen::Map<V>(source.at(ent.first).data(),
                                                                       source.at(ent.first).size())));
    return sink;
  }

}

#endif //HYBRIDPNM_EIGENSETGET_H
