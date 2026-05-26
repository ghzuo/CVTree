/*
 * Copyright (c) 2020  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-08 16:13:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2026-05-20 Wednesday 20:46:39
 */

#ifndef KIT_H
#define KIT_H

#include "fileOpt.h"
#include "info.h"
#include "jsonctl.h"
#include "ompOpt.h"
#include "stringOpt.h"

template <typename Key, typename Value, typename Hash = std::hash<Key>>
std::vector<std::pair<Key, Value>>
extractAndSort(std::unordered_map<Key, Value, Hash> &hashMap) {
  std::vector<std::pair<Key, Value>> result;
  result.reserve(hashMap.size());

  for (auto it = hashMap.begin(); it != hashMap.end();) {
    auto node = hashMap.extract(it++);
    result.emplace_back(std::move(node.key()), std::move(node.mapped()));
  }

  std::sort(result.begin(), result.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  return result;
};

// temporary index struct
template <typename T = size_t> struct Index {
  T value;
  Index(T &idx) : value(idx){};
  Index(size_t idx) : value(idx){};
  Index() : value(sentinel()){};
  bool isValid() const { return value != sentinel(); }
  static T invalid() { return sentinel(); }
  operator T() const { return value; }
  operator size_t() const { return static_cast<size_t>(value); }
  bool operator==(const Index &other) const { return value == other.value; }
  bool operator!=(const Index &other) const { return value != other.value; }

  friend ostream &operator<<(ostream &os, const Index &idx) {
    return os << idx.value;
  }
  friend istream &operator>>(istream &is, Index &idx) {
    return is >> idx.value;
  }

private:
  static constexpr T sentinel() { return numeric_limits<T>::max(); }
};
#endif // !KIT_H
