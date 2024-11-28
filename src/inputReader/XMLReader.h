#include <string>
#include "Simulation.h"


class XMLReader {
public:
    std::unique_ptr<Simulation> readXML(std::string fileName);

};

