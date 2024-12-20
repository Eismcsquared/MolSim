#include <string>
#include "Simulation.h"

/**
 * @brief The class provides the functionality to parse xml input files.
 */
class XMLReader {
public:
    /**
     * Read an xml input file.
     * @param particles The vector in which the particles should be stored.
     * @param fileName The name of the input file.
     * @return A Simulation object that contains all parameters and objects specified in the input file for the simulation.
     */
    static std::unique_ptr<Simulation> readXML(std::vector<Particle> &particles, std::string fileName);

};

