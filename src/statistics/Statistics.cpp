#include "Statistics.h"
#include "utils/ArrayUtils.h"
#include <utility>
#include <iomanip>
#include <fstream>
#include <spdlog/spdlog.h>

Statistics::Statistics(std::string file, int period, double from, double to, int numberBins, Axis profileAxis, Axis velocityAxis):
file(std::move(file)), period(period), from(from), to(to), numberBins(numberBins), profileAxis(profileAxis), velocityAxis(velocityAxis){}

std::string Statistics::getFile() const {
    return file;
}

int Statistics::getPeriod() const {
    return period;
}

double Statistics::getFrom() const {
    return from;
}

double Statistics::getTo() const {
    return to;
}

int Statistics::getNumberBins() const {
    return numberBins;
}

Axis Statistics::getProfileAxis() const {
    return profileAxis;
}

Axis Statistics::getVelocityAxis() const {
    return velocityAxis;
}

std::vector<double>
Statistics::aggregate(std::vector<Particle> &particles, const std::function<double(Particle&)>& f) {
    std::vector<double> result(numberBins, 0);
    double binSize = (to - from) / numberBins;
    for (Particle &particle : particles) {
        if (!particle.isStationary() && particle.isInDomain()) {
            int bin = std::floor((particle.getX()[profileAxis] - from) / binSize);
            if (bin >= 0 && bin < numberBins) {
                result[bin] += f(particle);
            }
        }
    }
    return result;
}

std::vector<double> Statistics::count(std::vector<Particle> &particles) {
    return aggregate(particles, [](Particle &particle) {return 1;});
}

std::vector<double> Statistics::density(std::vector<Particle> &particles) {
    return numberBins / (to - from) * count(particles);
}

std::vector<double> Statistics::velocity(std::vector<Particle> &particles) {
    const int axis = velocityAxis;
    std::vector<double> totalVelocity = aggregate(particles, [&axis](Particle &particle) {return particle.getV()[axis];});
    std::vector<double> number = count(particles);
    std::vector<double> result(numberBins, 0);
    for (int i = 0; i < numberBins; ++i) {
        if (number[i] >= 1) {
            result[i] = totalVelocity[i] / number[i];
        }
    }
    return result;
}

void Statistics::saveStatistics(std::vector<Particle> &particles, int iteration) {
    std::vector<double> rho = density(particles);
    std::vector<double> v = velocity(particles);
    std::stringstream fileName;
    fileName << file << "_" << std::setfill('0') << std::setw(4) << iteration << ".csv";
    std::ofstream outFile(fileName.str().c_str());
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (int i = 0; i < numberBins; ++i) {
        outFile << (i + 0.5) * (to - from) / numberBins + from << ", " << rho[i] << ", " << v[i] << std::endl;
    }
    outFile.close();
    spdlog::trace("The statistics of the system is successfully saved to {}", fileName.str());
}


