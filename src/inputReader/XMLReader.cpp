#include "XMLReader.h"
#include "inputReader/InputData.h"

std::unique_ptr<Simulation> XMLReader::readXML(std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        spdlog::error("Error: could not open file {}", fileName);
        exit(-1);
    }

    try {
        std::unique_ptr<InputData> input(simulation(file));
        std::unique_ptr<std::vector<Particle>> particles_ptr = std::make_unique<std::vector<Particle>>();
        std::vector<Particle>& particles = *particles_ptr;
    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
