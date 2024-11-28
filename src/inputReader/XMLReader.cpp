#include <fstream>
#include <memory>
#include <vector>
#include <spdlog/spdlog.h>
#include "XMLReader.h"
#include "inputReader/InputData.h"
#include "container/LinkedCellContainer.h"
#include "force/Force.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"
#include "container/DirectSumContainer.h"

std::unique_ptr<Simulation> XMLReader::readXML(std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        spdlog::error("Error: could not open file {}", fileName);
        exit(-1);
    }

    try {
        std::unique_ptr<InputData> input(simulation(file));
        std::unique_ptr<std::vector<Particle>> particles = std::make_unique<std::vector<Particle>>();
        for (auto p: input->objects().particles()) {
            particles->emplace_back(
                        std::array<double, 3>{p.position().x(), p.position().y(), p.position().z()},
                        std::array<double, 3>{p.velocity().x(), p.velocity().y(), p.velocity().z()},
                        p.mass()
                    );
        }
        for (auto c : input->objects().cuboids()) {
            Cuboid cuboid(
                        std::array<double, 3>{c.position().x(), c.position().y(), c.position().z()},
                        std::array<double, 3>{c.velocity().x(), c.velocity().y(), c.velocity().z()},
                        std::array<unsigned int, 3>{(unsigned int) c.size().x(), (unsigned int) c.size().y(), (unsigned int) c.size().z()},
                        c.mass(),
                        c.distance(),
                        c.brownVelocity(),
                        (int) c.brownDimension()
                    );
            cuboid.createParticles(*particles);
        }
        std::unique_ptr<Force> force;
        switch (input->parameters().force()) {
            case ForceType::value::gravitation:
                force = std::make_unique<GravitationalForce>();
                break;
            case ForceType::Lennard_Jones:
                force = std::make_unique<LennardJonesForce>();
                break;
        }
        std::unique_ptr<ParticleContainer> container;
        if (input->parameters().linked_cell().present()) {
            PositiveDoubleVector3 domainSize = input->parameters().linked_cell()->domain_size();
            PositiveDouble cutoffRadius = input->parameters().linked_cell()->cutoff_radius();
            container = std::make_unique<LinkedCellContainer>(
                    particles,
                    force,
                    std::array<double, 3>{domainSize.x(), domainSize.y(), domainSize.z()},
                    cutoffRadius);
        } else {
            container = std::make_unique<DirectSumContainer>(particles, force);
        }
        return std::make_unique<Simulation>(
                    container,
                    input->parameters().end_time(),
                    input->parameters().delta_t(),
                    input->parameters().output(),
                    input->parameters().format(),
                    input->parameters().frequency()
                );

    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
