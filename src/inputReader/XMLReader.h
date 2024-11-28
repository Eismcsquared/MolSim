#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <spdlog/spdlog.h>
#include "Simulation.h"

class XMLReader {
public:
    std::unique_ptr<Simulation> readXML(std::string fileName);

};

