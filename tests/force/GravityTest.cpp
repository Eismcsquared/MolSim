#include <gtest/gtest.h>
#include <vector>
#include <functional>
#include "spdlog/spdlog.h"
#include "container/DirectSumContainer.h"
#include "container/LinkedCellContainer.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"

class GravityTest: public ::testing::Test {
protected:

    // Both types of containers should be tested.
    std::vector<Particle> particles1;
    std::vector<Particle> particles2;
    std::unique_ptr<ParticleContainer> directSum;
    std::unique_ptr<ParticleContainer> linkedCell;
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

    void SetUp() override {

        std::unique_ptr<Force> g = std::make_unique<GravitationalForce>();
        std::unique_ptr<Force> l = std::make_unique<LennardJonesForce>();
        directSum = std::make_unique<DirectSumContainer>(particles1, g);
        linkedCell = std::make_unique<LinkedCellContainer>(particles2, l, std::array<double, 3>{15, 15, 1}, 3,
                                                           std::array<BoundaryCondition, 6>{REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING});
        directSum->setG(-10);
        linkedCell->setG(-10);
    }
};

// Test whether the motion of a single particle in the gravitation field corresponds to a parable.
TEST_F(GravityTest, OneParticle) {
    test_logger->info("GravityTest - One particle test");

    directSum->addParticle(Particle({1, 1, 0.5}, {2, 3, 0}, 1));
    directSum->simulate(0.5, 1e-5, "", "vtu", 10, false);

    EXPECT_LE(ArrayUtils::L2Norm(directSum->getParticles()[0].getX() - std::array<double, 3>{2, 1.25, 0.5}), 1e-5);
    EXPECT_LE(ArrayUtils::L2Norm(directSum->getParticles()[0].getV() - std::array<double, 3>{2, -2, 0}), 1e-5);

    linkedCell->addParticle(Particle({1, 1, 0.5}, {2, 3, 0}, 1));
    linkedCell->simulate(0.5, 1e-5, "", "vtu", 10, false);

    EXPECT_LE(ArrayUtils::L2Norm(linkedCell->getParticles()[0].getX() - std::array<double, 3>{2, 1.25, 0.5}), 1e-5);
    EXPECT_LE(ArrayUtils::L2Norm(linkedCell->getParticles()[0].getV() - std::array<double, 3>{2, -2, 0}), 1e-5);

    if (::testing::Test::HasFailure()) {
        test_logger->info("GravityTest - One particle test failed\n\n");
    } else {
        test_logger->info("GravityTest - One particle test passed\n\n");
    }
}

// Test the gravity in presence of interaction between particles. Since the exact motion is hard to predict, the test is based on energy conservation.
TEST_F(GravityTest, ManyBody) {
    test_logger->info("GravityTest - Many body test");

    auto potentialEnergy = [](std::vector<Particle> &particles, std::function<double(Particle &p1, Particle &p2)> &potential) {
        double E_pot = 0;
        for (int i = 0; i < particles.size(); ++i) {
            for (int j = i + 1; j < particles.size(); ++j) {
                E_pot += potential(particles[i], particles[j]);
            }
        }
        return E_pot;
    };

    auto kineticEnergy = [] (std::vector<Particle> &particles) {
        double E_kin = 0;
        for (Particle &p: particles) {
            E_kin += p.getM() * pow(ArrayUtils::L2Norm(p.getV()), 2) / 2;
        }
        return E_kin;
    };

    auto gravitationalPotential = [](std::vector<Particle> &particles, double g) {
        double E_g = 0;
        for (Particle &p: particles) {
            E_g += - p.getM() * g * p.getX()[1];
        }
        return E_g;
    };

    std::function<double(Particle &p1, Particle &p2)> gravitation = [](Particle &p1, Particle &p2) {
        return -p1.getM() * p2.getM() / ArrayUtils::L2Norm(p1.getX() - p2.getX());
    };

    std::function<double(Particle &p1, Particle &p2)> LennardJonesCutoff = [](Particle &p1, Particle &p2) {
        double r = std::min(ArrayUtils::L2Norm(p1.getX() - p2.getX()), 3.0);
        double a = pow(1 / r, 6);
        return 20 * (pow(a, 2) - a);
    };

    directSum->addParticle(Particle(std::array<double, 3>{7.5, 7.5, 0.5}, std::array<double, 3>{0, 0, 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{10, 7.5, 0.5}, std::array<double, 3>{0.2, 0, 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{8.75, 7.5 + 2.5 * sqrt(3) / 2, 0.5}, std::array<double, 3>{0.1, 0.1 * sqrt(3), 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{6.25, 7.5 + 2.5 * sqrt(3) / 2, 0.5}, std::array<double, 3>{-0.1, 0.1 * sqrt(3), 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{5, 7.5, 0.5}, std::array<double, 3>{-0.2, 0, 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{6.25, 7.5 - 2.5 * sqrt(3) / 2, 0.5}, std::array<double, 3>{-0.1, -0.1 * sqrt(3), 0}, 1));
    directSum->addParticle(Particle(std::array<double, 3>{8.75, 7.5 - 2.5 * sqrt(3) / 2, 0.5}, std::array<double, 3>{0.1, -0.1 * sqrt(3), 0}, 1));

    double energyBefore = potentialEnergy(directSum->getParticles(), gravitation) + kineticEnergy(directSum->getParticles()) +
            gravitationalPotential(directSum->getParticles(), -10);

    directSum->simulate(1, 1e-4, "", "vtu", 10, false);

    double energyAfter = potentialEnergy(directSum->getParticles(), gravitation) + kineticEnergy(directSum->getParticles()) +
                          gravitationalPotential(directSum->getParticles(), -10);
    EXPECT_NEAR(energyBefore, energyAfter, 1e-5);

    linkedCell->addParticle(Particle(std::array<double, 3>{7.5, 7.5, 0.5}, std::array<double, 3>{0, 0, 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{8.5, 7.5, 0.5}, std::array<double, 3>{0.2, 0, 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{8, 7.5 + sqrt(3) / 2, 0.5}, std::array<double, 3>{0.1, 0.1 * sqrt(3), 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{7, 7.5 + sqrt(3) / 2, 0.5}, std::array<double, 3>{-0.1, 0.1 * sqrt(3), 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{6.5, 7.5, 0.5}, std::array<double, 3>{-0.2, 0, 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{7, 7.5 - sqrt(3) / 2, 0.5}, std::array<double, 3>{-0.1, -0.1 * sqrt(3), 0}, 1));
    linkedCell->addParticle(Particle(std::array<double, 3>{8, 7.5 - sqrt(3) / 2, 0.5}, std::array<double, 3>{0.1, -0.1 * sqrt(3), 0}, 1));

    energyBefore = potentialEnergy(linkedCell->getParticles(), LennardJonesCutoff) + kineticEnergy(linkedCell->getParticles()) +
                          gravitationalPotential(linkedCell->getParticles(), -10);

    linkedCell->simulate(1, 1e-5, "", "vtu", 10, false);

    energyAfter = potentialEnergy(linkedCell->getParticles(), LennardJonesCutoff) + kineticEnergy(linkedCell->getParticles()) +
                         gravitationalPotential(linkedCell->getParticles(), -10);
    EXPECT_NEAR(energyBefore, energyAfter, 1e-5);

    if (::testing::Test::HasFailure()) {
        test_logger->info("GravityTest - Many body test failed\n\n");
    } else {
        test_logger->info("GravityTest - Many body test passed\n\n");
    }
}