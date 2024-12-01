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
            double r_z = p.position().z().present() ? p.position().z().get() : DoubleVector3::z_default_value();
            double v_z = p.velocity().z().present() ? p.velocity().z().get() : DoubleVector3::z_default_value();
            particles->emplace_back(
                        std::array<double, 3>{p.position().x(), p.position().y(), r_z},
                        std::array<double, 3>{p.velocity().x(), p.velocity().y(), v_z},
                        p.mass()
                    );
        }
        for (auto c : input->objects().cuboid()) {
            double r_z = c.position().z().present() ? c.position().z().get() : DoubleVector3::z_default_value();
            double v_z = c.velocity().z().present() ? c.velocity().z().get() : DoubleVector3::z_default_value();
            unsigned int n_z = c.size().z().present() ? c.size().z().get() : PositiveIntVector3::z_default_value();
            int dim = c.brownDimension().present() ? (int) c.brownDimension().get() : (int) CuboidType::brownDimension_default_value();
            Cuboid cuboid(
                        std::array<double, 3>{c.position().x(), c.position().y(), r_z},
                        std::array<double, 3>{c.velocity().x(), c.velocity().y(), v_z},
                        std::array<unsigned int, 3>{(unsigned int) c.size().x(), (unsigned int) c.size().y(), n_z},
                        c.mass(),
                        c.distance(),
                        c.brownVelocity(),
                        dim
                    );
            cuboid.createParticles(*particles);
        }
        std::unique_ptr<Force> force;
        if (input->parameters().gravitation().present()) {
            double g = input->parameters().gravitation().get().g().present() ? input->parameters().gravitation().get().g().get() : GravitationType::g_default_value();
            force = std::make_unique<GravitationalForce>(g);
        } else if (input->parameters().Lennard_Jones().present()) {
            double epsilon = input->parameters().Lennard_Jones().get().epsilon().present() ? input->parameters().Lennard_Jones().get().epsilon().get() : LennardJonesType::epsilon_default_value();
            double sigma = input->parameters().Lennard_Jones().get().sigma().present() ? input->parameters().Lennard_Jones().get().sigma().get() : LennardJonesType::sigma_default_value();
            force = std::make_unique<LennardJonesForce>(epsilon, sigma);
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
            BoundaryConditionType left = b.left().present() ? b.left().get() : BoundaryCondition3::left_default_value();
            BoundaryConditionType right = b.right().present() ? b.right().get() : BoundaryCondition3::right_default_value();
            BoundaryConditionType down = b.down().present() ? b.down().get() : BoundaryCondition3::down_default_value();
            BoundaryConditionType up = b.up().present() ? b.up().get() : BoundaryCondition3::up_default_value();
            BoundaryConditionType back = b.back().present() ? b.back().get() : BoundaryCondition3::back_default_value();
            BoundaryConditionType front = b.front().present() ? b.front().get() : BoundaryCondition3::front_default_value();
            std::array<BoundaryCondition, 6> boundaryCondition = {
                    boundaryConditionMap(left),
                    boundaryConditionMap(right),
                    boundaryConditionMap(down),
                    boundaryConditionMap(up),
                    boundaryConditionMap(back),
                    boundaryConditionMap(front)
            };
            double size_z = domainSize.z().present() ? domainSize.z().get() : PositiveDoubleVector3::z_default_value();
            container = std::make_unique<LinkedCellContainer>(
                    particles,
                    force,
                    std::array<double, 3>{domainSize.x(), domainSize.y(), size_z},
                    cutoffRadius,
                    boundaryCondition
            );
        } else {
            container = std::make_unique<DirectSumContainer>(particles, force);
        }
        double end_time = input->parameters().end_time().present() ? input->parameters().end_time().get() : SimulationParameters::end_time_default_value();
        double delta_t = input->parameters().delta_t().present() ? input->parameters().delta_t().get() : SimulationParameters::delta_t_default_value();
        std::string output = input->parameters().output().present() ? input->parameters().output().get() : SimulationParameters::output_default_value();
        std::string format = input->parameters().format().present() ? input->parameters().format().get() : SimulationParameters::format_default_value();
        unsigned int frequency = input->parameters().frequency().present() ? input->parameters().frequency().get() : SimulationParameters::frequency_default_value();

        return std::make_unique<Simulation>(
                    container,
                    end_time,
                    delta_t,
                    output,
                    format,
                    frequency
                );

    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
