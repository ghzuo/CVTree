#ifndef ARGPARSE_FACTORY_H
#define ARGPARSE_FACTORY_H

#include <argparse/argparse.hpp>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

// Define Parser as alias for convenience
// using Parser = argparse::ArgumentParser;
using json = nlohmann::json;

// Extended ArgumentParser class that provides access to argument map
class Parser : public argparse::ArgumentParser {
public:
  // Inherit constructors
  using argparse::ArgumentParser::ArgumentParser;

  // Get all argument names from the parser
  std::vector<std::string> get_long_arguments() const {
    std::vector<std::string> arguments;
    for (const auto &pair : this->m_argument_map) {
      if (pair.first.substr(0, 2) == "--") {
        arguments.emplace_back(pair.first);
      }
    }
    return arguments;
  };
};

// Parameter item structure to store argument information
struct ParamItem {
  std::string shortName;            // Short parameter name (e.g., "-i")
  std::string longName;             // Long parameter name (e.g., "--input")
  std::string helpText;             // Help text
  int nargs;                        // Number of arguments
  std::vector<std::string> choices; // Choices for the argument
  std::string jsonKey;              // Key name for JSON output

  // Default constructor
  ParamItem() : nargs(1) {}

  // Constructor
  ParamItem(const std::string &shortName, const std::string &longName,
            const std::string &helpText, int nargs = 1,
            const std::vector<std::string> &choices = {})
      : shortName(shortName), longName(longName), helpText(helpText),
        nargs(nargs), choices(choices),
        jsonKey(longName.substr(0, 2) == "--" ? longName.substr(2) : longName) {
  }
};

// Factory class for argparse management
class ArgparseFactory {
private:
  std::map<std::string, ParamItem>
      paramMap; // Parameter storage with jsonKey as key

public:
  // Add a parameter to the factory
  void addParam(const ParamItem &param) {
    if (paramMap.find(param.jsonKey) != paramMap.end()) {
      throw std::invalid_argument("Parameter with jsonKey '" + param.jsonKey +
                                  "' already exists");
    }
    paramMap[param.jsonKey] = param;
  }

  // Add parameters to an existing ArgumentParser
  void addParamsToParser(Parser &parser,
                         const std::vector<std::string> &jsonKeys) {
    // Add requested parameters to the parser
    for (const auto &jsonKey : jsonKeys) {
      auto it = paramMap.find(jsonKey);
      if (it == paramMap.end()) {
        throw std::invalid_argument("Parameter with jsonKey '" + jsonKey +
                                    "' not found in paramMap");
      }
      const ParamItem &param = it->second;
      auto &arg = parser.add_argument(param.shortName, param.longName)
                      .nargs(param.nargs);
      // Add choices if any
      if (!param.choices.empty()) {
        for (const auto &choice : param.choices) {
          arg.add_choice(choice);
        }
      }
      // add the help text
      arg.help(param.helpText);
    }
  }

  // Serialize ArgumentParser to JSON string
  std::string toJson(const Parser &parser) {
    json j = json::object();
    // Iterate through all parameters in the map
    const auto &arguments = parser.get_long_arguments();
    for (const auto &arg : arguments) {
      if (parser.is_used(arg)) {
        std::string jsonKey(arg.substr(2));
        auto it = paramMap.find(jsonKey);
        if (it->second.nargs == 0) {
          j[jsonKey] = true;
        } else {
          j[jsonKey] = parser.get<std::string>(arg);
        }
      }
    }
    return j.dump();
  }
};

#endif // ARGPARSE_FACTORY_H