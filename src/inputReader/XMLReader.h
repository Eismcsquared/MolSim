#include <string>
#include "Simulation.h"


class XMLReader {
public:
    static std::unique_ptr<Simulation> readXML(std::string fileName);

};

