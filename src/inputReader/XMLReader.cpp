#include <fstream>
#include <memory>
#include <vector>
#include <spdlog/spdlog.h>
#include <iostream>
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
        std::unique_ptr<InputData> input(simulation(file, xsd::cxx::tree::flags::dont_validate));
        std::unique_ptr<std::vector<Particle>> particles = std::make_unique<std::vector<Particle>>();
        for (auto p: input->objects().particle()) {
            particles->emplace_back(
                        std::array<double, 3>{p.position().x(), p.position().y(), p.position().z().get()},
                        std::array<double, 3>{p.velocity().x(), p.velocity().y(), p.velocity().z().get()},
                        p.mass()
                    );
        }
        for (auto c : input->objects().cuboid()) {
            Cuboid cuboid(
                        std::array<double, 3>{c.position().x(), c.position().y(), c.position().z().get()},
                        std::array<double, 3>{c.velocity().x(), c.velocity().y(), c.velocity().z().get()},
                        std::array<unsigned int, 3>{(unsigned int) c.size().x(), (unsigned int) c.size().y(), (unsigned int) c.size().z().get()},
                        c.mass(),
                        c.distance(),
                        c.brownVelocity(),
                        (int) c.brownDimension().get()
                    );
            cuboid.createParticles(*particles);
        }
        std::unique_ptr<Force> force;
        switch (input->parameters().force().get()) {
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
            BoundaryCondition3 b = input->parameters().linked_cell()->boundary_condition();
            auto boundaryConditionMap = [](const BoundaryConditionType& t) {
                switch (t) {
                    case BoundaryConditionType::value::reflection:
                        return REFLECTING;
                    case BoundaryConditionType::value::outflow:
                        return OUTFLOW;
                    case BoundaryConditionType::value::periodic:
                        return PERIODIC;
                }
            };
            std::array<BoundaryCondition, 6> boundaryCondition = {
                    boundaryConditionMap(b.left().get()),
                    boundaryConditionMap(b.right().get()),
                    boundaryConditionMap(b.down().get()),
                    boundaryConditionMap(b.up().get()),
                    boundaryConditionMap(b.back().get()),
                    boundaryConditionMap(b.front().get())
            };
            container = std::make_unique<LinkedCellContainer>(
                    particles,
                    force,
                    std::array<double, 3>{domainSize.x(), domainSize.y(), domainSize.z().get()},
                    cutoffRadius,
                    boundaryCondition
            );
        } else {
            container = std::make_unique<DirectSumContainer>(particles, force);
        }

        return std::make_unique<Simulation>(
                    container,
                    input->parameters().end_time().get(),
                    input->parameters().delta_t().get(),
                    input->parameters().output().get(),
                    input->parameters().format().get(),
                    input->parameters().frequency().get()
                );

    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
