#include <spdlog/spdlog.h>
#include <iomanip>
#include "StateWriter.h"


void StateWriter::saveState(std::vector<Particle> &particles, std::string fileName) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        spdlog::error("Error: Could not file {}", fileName);
        std::exit(-1);
    }

    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            for (double xi: p.getX()) {
                outFile << xi << " ";
            }
            for (double vi: p.getV()) {
                outFile << vi << " ";
            }
            outFile << p.getM() << " ";
            outFile << p.getType() << " ";
            outFile << p.getEpsilon() << " ";
            outFile << p.getSigma() << std::endl;
        }
    }

    outFile.close();
    spdlog::info("The state of the system is successfully saved to {}", fileName);
}
