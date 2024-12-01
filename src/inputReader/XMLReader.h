#include <string>
#include "Simulation.h"

/**
 * @brief The class provides the functionality to parse xml input files.
 */
class XMLReader {
public:
    static std::unique_ptr<Simulation> readXML(std::string fileName);

};

