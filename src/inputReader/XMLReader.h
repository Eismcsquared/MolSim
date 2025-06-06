#include <string>
#include "Simulation.h"
#include "InputData.h"

/**
 * @brief Represents data contains in a particle for parsing.
 */
struct ParticleData {
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    double mass;
    int type;
    double epsilon;
    double sigma;
};


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

    /**
     * Helper function that parse common data of all objects (Particle, Cuboid, Sphere, Membrane, Wall).
     * @param input Object that is read from the xml input file.
     * @return A struct that contains the data: position, velocity, mass, type, epsilon, sigma.
     */
    static ParticleData parseParticle(ParticleType& input);

};


