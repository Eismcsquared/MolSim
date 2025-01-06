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
    test_logger->info("ExternalForce - Force calculation test");

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
        test_logger->info("ExternalForce - Force calculation test failed\n\n");
    } else {
        test_logger->info("ExternalForce - Force calculation test passed\n\n");
    }
}
