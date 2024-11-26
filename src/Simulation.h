#include <memory>
#include "container/Container.h"

/**
 * @brief The class that manages a container of particles and stores the entire parameters for a simulation.
 */
class Simulation {
    /**
     * The objects in the simulation
     */
    std::unique_ptr<Container> container;
    /**
     * The end time of the simulation.
     */
    double endTime;
    /**
     * The time step of the simulation.
     */
    double deltaT;
    /**
     * The output format of the simulation data, either "vtu" or "xyz".
     */
    std::string outputFormat;
    /**
     * The name of the output file.
     */
    std::string outputFile;
    /**
     * The frequency of the output.
     */
    unsigned int outputFrequency;
    /**
     * Represents whether output should be activated, useful for benchmarking.
     */
    bool output;
    /**
     * represents whether the Newton's third law should be applied in the force calculations.
     */
    bool newton3;

};

