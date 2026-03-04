#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <boost/functional/hash.hpp>

/*
Hash function for std::pair

Syntax taken from: https://lists.isocpp.org/std-discussion/2020/12/0937.php
Using boost's hash_combine:
https://www.boost.org/doc/libs/1_51_0/doc/html/hash/combine.html

*/

namespace std {

template <typename T1, typename T2>
struct hash<pair<T1, T2>> {
  size_t operator()(const pair<T1, T2> &input_pair) const {
    size_t seed = 0;
    boost::hash_combine(seed, input_pair.first);
    boost::hash_combine(seed, input_pair.second);
    return seed;
    // return std::hash<T1>()(input_pair.first) ^
    // std::hash<T2>()(input_pair.second);
  }
};
}  // namespace std

/**
 * @brief Search file paths recursively.
 *
 * @details If a directory path is given, it returns all the file paths in the
 * directory and subdirectories. If a file path is given, it returns the file
 * path.
 * @param path Directory or file path.
 * @return Returns a vector of found file paths.
 */
inline std::vector<std::filesystem::path> find_file_paths(
    const std::filesystem::path path) {
  std::vector<std::filesystem::path> paths;
  if (std::filesystem::is_regular_file(std::filesystem::path(path))) {
    paths.emplace_back(path);
  } else {
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(path)) {
      if (entry.is_regular_file()) paths.emplace_back(entry.path());
    }
  }
  return paths;
}