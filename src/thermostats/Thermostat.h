#include <vector>
#include "body/Particle.h"

/**
 * @brief This class represent a thermostat that can adjust the temperature of a system.
 */
class Thermostat {

private:
    /**
     * The target temperature.
     */
    double target_T;

    /**
     * The application period of the thermostat in terms of iteration number.
     */
    int periode;

    /**
     * The maximal change of the temperature in one step.
     */
    double maxDelta;

    /**
     * The dimension of the system.
     */
     int dimension;

public:

    /**
     * Constructor.
     * @param targetT The target temperature.
     * @param periode The application period of the thermostat.
     * @param maxDelta The maximal change of temperature in one step.
     * @param dimension The dimension of the simulation.
     */
    Thermostat(double targetT, int periode, double maxDelta, int dimension);

    /**
     * Apply the thermostat.
     * @param particles The particles which the thermostat should be applied to.
     */
    void apply(std::vector<Particle> &particles) const;

    /**
     * Getter for the application periode.
     * @return The application periode.
     */
    int getPeriode() const;
};

