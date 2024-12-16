#include <fstream>
#include <memory>
#include <vector>
#include <spdlog/spdlog.h>
#include <limits>
#include <iostream>
#include "XMLReader.h"
#include "inputReader/InputData.h"
#include "container/LinkedCellContainer.h"
#include "force/Force.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"
#include "container/DirectSumContainer.h"
#include "body/Sphere.h"
#include "inputReader/StateReader.h"

std::unique_ptr<Simulation> XMLReader::readXML(std::vector<Particle> &particles, std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        spdlog::error("Error: could not open file {}", fileName);
        exit(-1);
    }

    try {
        std::unique_ptr<InputData> input(simulation(file, xsd::cxx::tree::flags::dont_validate));



        for (auto p: input->objects().particle()) {
            double r_z = p.position().z().present() ? p.position().z().get() : PositiveDoubleVector3::z_default_value() / 2;
            double v_z = p.velocity().z().present() ? p.velocity().z().get() : DoubleVector3::z_default_value();
            int type = p.type().present() ? p.type().get() : ParticleType::type_default_value();
            double epsilon = p.epsilon().present() ? p.epsilon().get() : ParticleType::epsilon_default_value();
            double sigma = p.sigma().present() ? p.sigma().get() : ParticleType::sigma_default_value();
            particles.emplace_back(
                        std::array<double, 3>{p.position().x(), p.position().y(), r_z},
                        std::array<double, 3>{p.velocity().x(), p.velocity().y(), v_z},
                        p.mass(),
                        type,
                        epsilon,
                        sigma
                    );
        }

        int dimension = static_cast<int>(input->parameters().dimension().present() ? input->parameters().dimension().get() : InputData::parameters_type::dimension_default_value());

        for (auto c : input->objects().cuboid()) {
            double r_z = c.position().z().present() ? c.position().z().get() : PositiveDoubleVector3::z_default_value() / 2;
            double v_z = c.velocity().z().present() ? c.velocity().z().get() : DoubleVector3::z_default_value();
            unsigned int n_z = c.size().z().present() ? c.size().z().get() : PositiveIntVector3::z_default_value();
            int type = c.type().present() ? c.type().get() : ParticleType::type_default_value();
            double epsilon = c.epsilon().present() ? c.epsilon().get() : CuboidType::epsilon_default_value();
            double sigma = c.sigma().present() ? c.sigma().get() : CuboidType::sigma_default_value();

            double brown_vel;
            if (input->parameters().thermostat().present()) {
                brown_vel = sqrt(input->parameters().thermostat()->initial_T() / c.mass());
            } else if (c.brown_velocity().present()) {
                brown_vel = c.brown_velocity().get();
            } else {
                brown_vel = CuboidType::brown_velocity_default_value();
            }

            Cuboid cuboid(
                        std::array<double, 3>{c.position().x(), c.position().y(), r_z},
                        std::array<double, 3>{c.velocity().x(), c.velocity().y(), v_z},
                        std::array<unsigned int, 3>{(unsigned int) c.size().x(), (unsigned int) c.size().y(), n_z},
                        c.mass(),
                        c.distance(),
                        brown_vel,
                        dimension,
                        type,
                        epsilon,
                        sigma
                    );

            cuboid.createParticles(particles);

        }

        for (auto s : input->objects().sphere()) {
            double r_z = s.center().z().present() ? s.center().z().get() : PositiveDoubleVector3::z_default_value() / 2;
            double v_z = s.velocity().z().present() ? s.velocity().z().get() : DoubleVector3::z_default_value();
            int radius = static_cast<int>(s.radius());
            int type = s.type().present() ? s.type().get() : ParticleType::type_default_value();
            double epsilon = s.epsilon().present() ? s.epsilon().get() : SphereType::epsilon_default_value();
            double sigma = s.sigma().present() ? s.sigma().get() : SphereType::sigma_default_value();

            double brown_vel;
            if (input->parameters().thermostat().present()) {
                brown_vel = sqrt(input->parameters().thermostat()->initial_T() / s.mass());
            } else if (s.brown_velocity().present()) {
                brown_vel = s.brown_velocity().get();
            } else {
                brown_vel = SphereType::brown_velocity_default_value();
            }

            Sphere sphere(
                        std::array<double, 3>{s.center().x(), s.center().y(), r_z},
                        std::array<double, 3>{s.velocity().x(), s.velocity().y(), v_z},
                        radius,
                        s.mass(),
                        s.distance(),
                        brown_vel,
                        dimension,
                        type,
                        epsilon,
                        sigma
                    );
            sphere.createParticles(particles);
        }

        for (auto f:  input->objects().load()) {
            StateReader::loadState(particles, f);
        }

        std::unique_ptr<Force> force;
        ForceType f = input->parameters().force().present() ? input->parameters().force().get() : SimulationParameters::force_default_value();
        switch (f) {
            case ForceType::value::gravitation:
                force = std::make_unique<GravitationalForce>();
                break;
            case ForceType::value::Lennard_Jones:
                force = std::make_unique<LennardJonesForce>();
                break;
        }

        std::unique_ptr<ParticleContainer> container;
        if (input->parameters().linked_cell().present()) {
            PositiveDoubleVector3 domainSize = input->parameters().linked_cell()->domain_size();
            PositiveDouble cutoffRadius = input->parameters().linked_cell()->cutoff_radius();
            auto b = input->parameters().linked_cell()->boundary_condition();
            auto boundaryConditionMap = [](const BoundaryConditionType& t) {
                switch (t) {
                    case BoundaryConditionType::value::reflect:
                        return REFLECTING;
                    case BoundaryConditionType::value::outflow:
                        return OUTFLOW;
                    case BoundaryConditionType::value::periodic:
                        return PERIODIC;
                }
                return OUTFLOW;
            };
            BoundaryConditionType left = b.present() && b.get().left().present() ? b.get().left().get() : BoundaryCondition3::left_default_value();
            BoundaryConditionType right = b.present() && b.get().right().present() ? b.get().right().get() : BoundaryCondition3::right_default_value();
            BoundaryConditionType down = b.present() && b.get().down().present() ? b.get().down().get() : BoundaryCondition3::down_default_value();
            BoundaryConditionType up = b.present() && b.get().up().present() ? b.get().up().get() : BoundaryCondition3::up_default_value();
            BoundaryConditionType back = b.present() && b.get().back().present() ? b.get().back().get() : BoundaryCondition3::back_default_value();
            BoundaryConditionType front = b.present() && b.get().front().present() ? b.get().front().get() : BoundaryCondition3::front_default_value();
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

        if (input->parameters().g().present()) {
            container->setG(input->parameters().g().get());
        }

        if (input->parameters().thermostat().present()) {
            auto ts = input->parameters().thermostat().get();
            double target_T = ts.target_T().present() ? ts.target_T().get() : ts.initial_T();
            double maxChange = ts.max_delta().present() ? static_cast<double>(ts.max_delta().get()) : std::numeric_limits<double>::infinity();
            std::unique_ptr<Thermostat> thermostat = std::make_unique<Thermostat>(
                    target_T,
                    ts.periode(),
                    maxChange,
                    dimension
                );
            container->setThermostat(thermostat);
        }

        double end_time = input->parameters().end_time().present() ? input->parameters().end_time().get() : SimulationParameters::end_time_default_value();
        double delta_t = input->parameters().delta_t().present() ? input->parameters().delta_t().get() : SimulationParameters::delta_t_default_value();
        std::string output = input->parameters().output().present() ? input->parameters().output().get() : SimulationParameters::output_default_value();
        std::string format = input->parameters().format().present() ? input->parameters().format().get() : SimulationParameters::format_default_value();
        unsigned int frequency = input->parameters().frequency().present() ? input->parameters().frequency().get() : SimulationParameters::frequency_default_value();

        std::unique_ptr<Simulation> simulation =  std::make_unique<Simulation>(
                    container,
                    end_time,
                    delta_t,
                    output,
                    format,
                    frequency
                );

        if (input->parameters().store().present()) {
            simulation->setCheckpointing(input->parameters().store().get());
            simulation->setSaveOutput(false);
        }

        return simulation;

    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
