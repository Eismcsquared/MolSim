#include <functional>
#include "body/Particle.h"

/**
 * @brief This class computes and saves statistics for the nanoflow simulation: The density and velocity profile.
 */
enum Axis {
    X, Y, Z
};

class Statistics {
private:

    /**
     * The file name for statistics output.
     */
    std::string file;

    /**
     * The period (number of iterations) of statistics output.
     */
    int period;

    /**
     * The coordinate of the start position.
     */
    double from;

    /**
     * The coordinate of the end position.
     */
    double to;

    /**
     * The number of bins, which the range [from, to) is divided into.
     */
    int numberBins;

    /**
     * The axis along which the profile is computed.
     */
    Axis profileAxis;

    /**
     * The axis along which the velocity is considered.
     */
    Axis velocityAxis;

public:

    /**
     * Constructors.
     * @param file The output file.
     * @param period The period of output.
     * @param from The start coordinate.
     * @param to The end position.
     * @param numberBins The number ofd bins, which the range [from, to) is divided into.
     * @param profileAxis The axis for profiles. Default: X.
     * @param velocityAxis The velocity component of interest. Default: Y.
     */
    Statistics(std::string file, int period, double from, double to, int numberBins, Axis profileAxis = X, Axis velocityAxis = Y);

    /**
     * Getter for the period.
     * @return The output period in number of iterations.
     */
    double getPeriod() const;

    /**
     * Aggregate particles into bins according to a given aggregation function.
     * @param particles The particles.
     * @param f The aggregation function.
     * @return A vector containing the sum over aggregation function values of the particles in each bin.
     */
    std::vector<double> aggregate(std::vector<Particle> &particles, const std::function<double(Particle&)>& f);

    /**
     * Count particles in each bin.
     * @param particles The particles to count.
     * @return A vector containing the number of particles in each bin.
     */
    std::vector<double> count(std::vector<Particle> &particles);

    /**
     * Compute the density profile.
     * @param particles The particles from which the density should be computed.
     * @return A vector containing the average particle number per length unit in each bin.
     */
    std::vector<double> density(std::vector<Particle> &particles);

    /**
     * Compute the velocity profile.
     * @param particles The particles from which the density should be computed.
     * @return A vector containing the average velocity in each bin.
     */
    std::vector<double> velocity(std::vector<Particle> &particles);

    /**
     * Write statistics into a csv file. Each row contains the following information, separated by kommas.
     * The coordinate of the center of a bin, the average length density of the bin, the average velocity in the bin.
     * @param particles The particles that should be analyzed.
     * @param iteration The current iteration number.
     */
    void saveStatistics(std::vector<Particle> &particles, int iteration);

};

