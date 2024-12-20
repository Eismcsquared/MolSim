#include <fstream>
#include <sstream>
#include <spdlog/spdlog.h>
#include "StateReader.h"

void StateReader::loadState(std::vector<Particle> &particles, std::string fileName) {
    std::ifstream inFile(fileName);

    if (!inFile.is_open()) {
        spdlog::error("Error: Could not open file {}", fileName);
        std::exit(-1);
    }

    std::array<double, 3> x{};
    std::array<double, 3> v{};
    double m;
    int type;
    double epsilon;
    double sigma;

    std::string line;
    while (std::getline(inFile, line)) {
        spdlog::trace("Read line: {}", line);
        std::istringstream dataStream(line);
        for (double &xi: x) {
            dataStream >> xi;
        }
        for (double &vi: v) {
            dataStream >> vi;
        }
        dataStream >> m;
        dataStream >> type;
        dataStream >> epsilon;
        dataStream >> sigma;

        particles.emplace_back(x, v, m, type, epsilon, sigma);
    }
}
