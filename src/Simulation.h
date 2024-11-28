#include <memory>
#include "container/ParticleContainer.h"

/**
 * @brief The class that manages a container of particles and stores the entire parameters for a simulation.
 */
class Simulation {
    /**
     * The objects in the simulation
     */
    std::unique_ptr<ParticleContainer> container;
    /**
     * The end time of the simulation.
     */
    double endTime;
    /**
     * The time step of the simulation.
     */
    double deltaT;

     /**
     * The name of the output file.
     */
    std::string outputFile;

    /**
     * The output format of the simulation data, either "vtu" or "xyz".
     */
    std::string outputFormat;

    /**
     * The frequency of the output.
     */
    unsigned int outputFrequency;
    /**
     * Represents whether output should be activated, useful for benchmarking.
     */
    bool saveOutput;
    /**
     * represents whether the Newton's third law should be applied in the force calculations.
     */
    bool newton3;
public:

    /**
     * Constructor.
     * @param container The container that contains the particles for the simulation.
     * @param endTime The duration of the simulation.
     * @param deltaT The time step of the simulation.
     * @param outputFormat The format of the output file.
     * @param outputFile The name of the output file.
     * @param outputFrequency The frequency of output.
     */
    Simulation(std::unique_ptr<ParticleContainer> &container, double endTime, double deltaT,
               std::string outputFile, std::string outputFormat, unsigned int outputFrequency);

    /**
     * Setter for the end time.
     * @param endTime The new end time.
     */
    void setEndTime(double endTime);

    /**
     * Setter for the time step.
     * @param endTime The new time step.
     */
    void setDeltaT(double deltaT);

    /**
     * Setter for the format of the output.
     * @param endTime The new output format, either "vtu" or "xyz".
     */
    void setOutputFormat(const std::string &outputFormat);

    /**
     * Setter for the name of the output file.
     * @param outputFile The new output file name.
     */
    void setOutputFile(const std::string &outputFile);

    /**
     * Setter for the frequency of output
     * @param outputFrequency The new output frequency.
     */
    void setOutputFrequency(unsigned int outputFrequency);

    /**
     * Setter for the flag of saving output.
     * @param saveOutput The new flag.
     */
    void setSaveOutput(bool saveOutput);

    /**
     * Setter for the flag of applying the Newton's third law.
     * @param saveOutput The new flag.
     */
    void setNewton3(bool newton3);

     /**
      * Run the simulation.
      */
    void run();

};

