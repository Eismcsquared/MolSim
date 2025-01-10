#include <fstream>
#include <memory>
#include <vector>
#include <spdlog/spdlog.h>
#include <limits>
#include "XMLReader.h"
#include "container/LinkedCellContainer.h"
#include "force/Force.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"
#include "force/MembraneForce.h"
#include "container/DirectSumContainer.h"
#include "body/Sphere.h"
#include "body/Membrane.h"
#include "inputReader/StateReader.h"

ParticleData XMLReader::parseParticle(ParticleType &input) {
    std::array<double, 3> r{input.position().x(), input.position().y(), input.position().z().present() ? input.position().z().get() : PositiveDoubleVector3::z_default_value() / 2};
    std::array<double, 3> v{0, 0, 0};
    if (input.velocity().present()) {
        v[0] = input.velocity().get().x();
        v[1] = input.velocity().get().y();
        if (input.velocity().get().z().present()) {
            v[2] = input.velocity().get().z().get();
        }
    }
    double m = input.mass();
    int type = input.type().present() ? input.type().get() : ParticleType::type_default_value();
    double epsilon = input.epsilon().present() ? input.epsilon().get() : ParticleType::epsilon_default_value();
    double sigma = input.sigma().present() ? input.sigma().get() : ParticleType::sigma_default_value();
    return ParticleData{r, v, m, type, epsilon, sigma};
}

std::unique_ptr<Simulation> XMLReader::readXML(std::vector<Particle> &particles, std::string fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        spdlog::error("Error: could not open file {}", fileName);
        exit(-1);
    }

    try {
        std::unique_ptr<InputData> input(simulation(file, xsd::cxx::tree::flags::dont_validate));

        for (auto checkpointingFile:  input->objects().load()) {
            StateReader::loadState(particles, checkpointingFile);
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
            case ForceType::value::membrane:
                force = std::make_unique<MembraneForce>();
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


        for (auto p: input->objects().particle()) {
            ParticleData data = parseParticle(p);
            container->addParticle(
                    Particle(
                        data.position,
                        data.velocity,
                        data.mass,
                        data.type,
                        data.epsilon,
                        data.sigma
                    )
            );
        }

        int dimension = static_cast<int>(input->parameters().dimension().present() ? input->parameters().dimension().get() : InputData::parameters_type::dimension_default_value());

        for (auto c : input->objects().cuboid()) {

            ParticleData data = parseParticle(c);
            unsigned int n_z = c.size().z().present() ? c.size().z().get() : PositiveIntVector3::z_default_value();
            double brown_vel;
            if (input->parameters().thermostat().present()) {
                brown_vel = sqrt(input->parameters().thermostat()->initial_T() / c.mass());
            } else if (c.brown_velocity().present()) {
                brown_vel = c.brown_velocity().get();
            } else {
                brown_vel = CuboidType::brown_velocity_default_value();
            }

            Cuboid cuboid(
                        data.position,
                        data.velocity,
                        std::array<unsigned int, 3>{(unsigned int) c.size().x(), (unsigned int) c.size().y(), n_z},
                        data.mass,
                        c.distance(),
                        brown_vel,
                        dimension,
                        data.type,
                        data.epsilon,
                        data.sigma
                    );

            container->addCluster(cuboid);

        }

        for (auto s : input->objects().sphere()) {

            ParticleData data = parseParticle(s);
            int radius = static_cast<int>(s.radius());
            double brown_vel;
            if (input->parameters().thermostat().present()) {
                brown_vel = sqrt(input->parameters().thermostat()->initial_T() / s.mass());
            } else if (s.brown_velocity().present()) {
                brown_vel = s.brown_velocity().get();
            } else {
                brown_vel = SphereType::brown_velocity_default_value();
            }

            Sphere sphere(
                        data.position,
                        data.velocity,
                        radius,
                        data.mass,
                        s.distance(),
                        brown_vel,
                        dimension,
                        data.type,
                        data.epsilon,
                        data.sigma
                    );
            container->addCluster(sphere);
        }

        for (auto m : input->objects().membrane()) {
            ParticleData data = parseParticle(m);
            double k = m.k().present() ? m.k().get() : MembraneType::k_default_value();
            double r0 = m.r0().present() ? m.r0().get() : MembraneType::r0_default_value();
            double brown_vel;
            if (input->parameters().thermostat().present()) {
                brown_vel = sqrt(input->parameters().thermostat()->initial_T() / m.mass());
            } else if (m.brown_velocity().present()) {
                brown_vel = m.brown_velocity().get();
            } else {
                brown_vel = CuboidType::brown_velocity_default_value();
            }

            Membrane membrane(
                    data.position,
                    data.velocity,
                    std::array<unsigned int, 2>{static_cast<unsigned int>(m.size().x()), static_cast<unsigned int>(m.size().y())},
                    data.mass,
                    m.distance(),
                    brown_vel,
                    dimension,
                    data.type,
                    data.epsilon,
                    data.sigma,
                    k,
                    r0
            );

            unsigned int oldSize = particles.size();

            container->addCluster(membrane);

            for (auto ef: m.external_force()) {
                double fx = ef.force().x().present() ? ef.force().x().get() : OptionalDoubleVector3::x_default_value();
                double fy = ef.force().y().present() ? ef.force().y().get() : OptionalDoubleVector3::y_default_value();
                double fz = ef.force().z().present() ? ef.force().z().get() : OptionalDoubleVector3::z_default_value();
                for (auto index: ef.index()) {
                    container->addExternalForce(oldSize + (index.x() - 1) * m.size().y() + (index.y() - 1), {fx, fy, fz}, ef.until());
                }
            }
        }

        for (auto w : input->objects().wall()) {

            ParticleData data = parseParticle(w);
            unsigned int n_z = w.size().z().present() ? w.size().z().get() : PositiveIntVector3::z_default_value();
            Cuboid wall(
                    data.position,
                    data.velocity,
                    std::array<unsigned int, 3>{static_cast<unsigned int>(w.size().x()), static_cast<unsigned int>(w.size().y()), n_z},
                    data.mass,
                    w.distance(),
                    0,
                    dimension,
                    data.type,
                    data.epsilon,
                    data.sigma,
                    true
            );
            container->addCluster(wall);
        }


        if (input->parameters().start_time().present()) {
            container->setT(input->parameters().start_time().get());
        }

        if (input->parameters().g().present()) {
            OptionalDoubleVector3 g = input->parameters().g().get();
            double gx = g.x().present() ? g.x().get() : OptionalDoubleVector3::x_default_value();
            double gy = g.y().present() ? g.y().get() : OptionalDoubleVector3::y_default_value();
            double gz = g.z().present() ? g.z().get() : OptionalDoubleVector3::z_default_value();
            container->setG({gx, gy, gz});
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

        if (output.empty()) {
            simulation->setSaveOutput(false);
        }

        if (input->parameters().store().present()) {
            simulation->setCheckpointingFile(input->parameters().store().get());
        }

        return simulation;

    } catch (const xml_schema::exception& e) {
        spdlog::error("XML parsing error: {}", e.what());
        std::exit(-1);
    }


}
