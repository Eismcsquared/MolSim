#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include "force/MembraneForce.h"
#include "container/LinkedCellContainer.h"
#include "body/Membrane.h"

class MembraneForceTest : public ::testing::Test {
protected:
    std::unique_ptr<Force> force;
    std::vector<Particle> particles;
    std::unique_ptr<ParticleContainer> pc;
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

    void SetUp() override {
        force = std::make_unique<MembraneForce>();
        pc = std::make_unique<LinkedCellContainer>(particles, force, std::array<double, 3>{10, 10, 10},3,
                                                   std::array<BoundaryCondition, 6>{OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW});
    }
};

// Test whether the membrane force is calculated for both neighboring and non-neighboring particles. For non-neighboring particles,
// cases distance < 2^(1/6)sigma and > 2^(1/6)sigma are considered.
TEST_F(MembraneForceTest, ForceCalculation) {
    test_logger->info("MembraneForce - Force calculation test");

    Membrane m({1, 1, 2}, {0, 0, 0}, {3, 3}, 1, 1, 0, 3, 0, 5, 1, 10, 1.5);
    // should interact with the central particle of the membrane.
    pc->addParticle(Particle({2, 2, 3.1}, {0, 0, 0}, 1));
    // no interaction with the membrane since too far.
    pc->addParticle(Particle({2, 2, 0.85}, {0, 0, 0}, 1));
    pc->addCluster(m);
    pc->updateF();

    EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[0].getF() - std::array<double, 3>{0, 0, 7.9404769491203128}), 1e-12);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[1].getF() - std::array<double, 3>{0, 0, 0}), 1e-12);
    for (int i = 2; i < 11; ++i) {
        if (ArrayUtils::L2Norm(pc->getParticles()[i].getX() - std::array<double, 3>{2, 2, 2}) < 1e-12) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -7.9404769491203128}), 1e-12);
        } else if (std::abs(ArrayUtils::L2Norm(pc->getParticles()[i].getX() - std::array<double, 3>{2, 2, 2}) - 1) < 1e-12) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - 15 * (pc->getParticles()[i].getX() - std::array<double, 3>{2, 2, 2})), 1e-12);
        } else {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - 10 * (pc->getParticles()[i].getX() - std::array<double, 3>{2, 2, 2})), 1e-12);
        }
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("MembraneForce - Force calculation test failed\n\n");
    } else {
        test_logger->info("MembraneForce - Force calculation test passed\n\n");
    }
}

// Test whether external forces are applied correctly to the corresponding particles until the specified time.
TEST_F(MembraneForceTest, ExternalForce) {
    test_logger->info("MembraneForce - External force test");

    Membrane m({1, 1, 1}, {0, 0, 0}, {3, 3}, 1, 1, 0, 3, 0, 5, 1, 10, 1);
    pc->addCluster(m);
    pc->addExternalForce(0, {-1, 1, 3}, 10);
    pc->addExternalForce(6, {0, 0, 2}, 15);
    pc->setG({0, 0, -10});
    pc->updateF();

    for (int i = 0; i < 9; ++i) {
        if (i == 0) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{-1, 1, -7}), 1e-12);
        } else if (i == 6) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -8}), 1e-12);
        } else {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -10}), 1e-12);
        }
    }

    // external forces on particle 0 and 6 apply.
    pc->setT(9.99);
    pc->updateF();

    for (int i = 0; i < 9; ++i) {
        if (i == 0) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{-1, 1, -7}), 1e-12);
        } else if (i == 6) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -8}), 1e-12);
        } else {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -10}), 1e-12);
        }
    }

    // external force on particle 0 no more applies.
    pc->setT(10.01);
    pc->updateF();

    for (int i = 0; i < 9; ++i) {
        if (i == 6) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -8}), 1e-12);
        } else {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -10}), 1e-12);
        }
    }

    // external force on particle 6 still applies.
    pc->setT(14.99);
    pc->updateF();

    for (int i = 0; i < 9; ++i) {
        if (i == 6) {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -8}), 1e-12);
        } else {
            EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -10}), 1e-12);
        }
    }

    // no more external force applies.
    pc->setT(15.01);
    pc->updateF();

    for (int i = 0; i < 9; ++i) {
        EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[i].getF() - std::array<double, 3>{0, 0, -10}), 1e-12);
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("MembraneForce - External force test failed\n\n");
    } else {
        test_logger->info("MembraneForce - External force test passed\n\n");
    }
}

// Compare the simulation of a 1D chain consisting of three particles with its analytical solution.
// With the following parameter:
// x_1(0) = 1.5, x_2(0) = 2.5, x_3(0) = 5, v_1(0) = 0, v_2(0) = 0, v_3(0) = 0, r_0 = 2, k = 100, m = 1
// The analytical solution is given as
// x_1(t) = 1 + 0.25 * cos(10 * t) + 0.25 * cos(10 * sqrt(3) * t)
// x_2(t) = 3 - 0.5 * cos(10 * sqrt(3) * t)
// x_3(t) = 5 - 0.25 * cos(10 * t) + 0.25 * cos(10 * sqrt(3) * t)
TEST_F(MembraneForceTest, Anaytical) {
    test_logger->info("MembraneForce - Analytical test");

    Particle p1({1.5, 1, 1}, {0, 0, 0}, 1);
    Particle p2({2.5, 1, 1}, {0, 0, 0}, 1);
    Particle p3({5, 1, 1}, {0, 0, 0}, 1);

    p1.addNeighbour(p2);
    p2.addNeighbour(p1);
    p2.addNeighbour(p3);
    p3.addNeighbour(p2);

    pc->addParticle(p1);
    pc->addParticle(p2);
    pc->addParticle(p3);

    for (Particle &p : pc->getParticles()) {
        p.setK(100);
        p.setR0(2);
    }

    auto solution1 = [](double t) { return 1 + 0.25 * std::cos(10 * t) + 0.25 * std::cos(10 * std::sqrt(3) * t); };
    auto solution2 = [](double t) { return 3 - 0.5 * std::cos(10 * std::sqrt(3) * t); };
    auto solution3 = [](double t) { return 5 - 0.25 * std::cos(10 * t) + 0.25 * std::cos(10 * std::sqrt(3) * t); };

    for (int i = 1; i <= 10; ++i) {
        pc->simulate(0.05 * i, 1e-5, "", "", 10, false);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[0].getX() - std::array<double, 3>{solution1(pc->getT()), 1, 1}), 1e-6);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[1].getX() - std::array<double, 3>{solution2(pc->getT()), 1, 1}), 1e-6);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(pc->getParticles()[2].getX() - std::array<double, 3>{solution3(pc->getT()), 1, 1}), 1e-6);
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("MembraneForce - Analytical test failed\n\n");
    } else {
        test_logger->info("MembraneForce - Analytical test passed\n\n");
    }
}