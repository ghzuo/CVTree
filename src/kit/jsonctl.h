/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-15 Monday 08:10:07
 */

#ifndef JSONCTL
#define JSONCTL

#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
using Json = nlohmann::json;
using namespace std;

struct JsonCtl {
  bool format{false};

  void operator()(const string &, const string &);
  void operator()(ostream &, const string &);
};

// read a value from json object
template <typename T>
T getJsonValue(const Json &jsonObj, const string &key, const T &defaultValue) {
  try {
    if (jsonObj.contains(key)) {
      return jsonObj[key].get<T>();
    } else {
      return defaultValue;
    }
  } catch (const Json::type_error &e) {
    cerr << "Error: Type conversion exception when reading JSON key '" << key
         << "': " << e.what() << endl;
    exit(1);
  } catch (const exception &e) {
    cerr << "Error: Exception (" << typeid(e).name()
         << ") when reading JSON key '" << key << "': " << e.what() << endl;
    exit(1);
  }
}

#endif // !JSONCTL

#ifndef THEJSON
#define THEJSON
extern JsonCtl theJson;
#endif //.!THEJSON
