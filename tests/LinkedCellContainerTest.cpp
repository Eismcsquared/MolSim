#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include "container/LinkedCellContainer.h"
#include "force/LennardJonesForce.h"
#include "Logger.h"

class LinkedCellContainerTest: public ::testing::Test {
protected:
    std::unique_ptr<LinkedCellContainer> container2D;
    std::unique_ptr<LinkedCellContainer> container3D;

    void SetUp() override {
        std::unique_ptr<std::vector<Particle>> p1 = std::make_unique<std::vector<Particle>>();
        std::unique_ptr<std::vector<Particle>> p2 = std::make_unique<std::vector<Particle>>();
        std::unique_ptr<Force> f1 = std::make_unique<LennardJonesForce>();
        std::unique_ptr<Force> f2 = std::make_unique<LennardJonesForce>();
        container2D = std::make_unique<LinkedCellContainer>(p1, f1, std::array<double, 3>{15, 15, 1}, 3,
                                          std::array<BoundaryCondition, 6>{OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW});
        container3D = std::make_unique<LinkedCellContainer>(p2, f2, std::array<double, 3>{100, 100, 50}, 3,
                                          std::array<BoundaryCondition, 6>{REFLECTING, REFLECTING, REFLECTING, REFLECTING, OUTFLOW, OUTFLOW});
        test_logger -> info("LinkedCellContainer created");
    }

    void TearDown() override {
        test_logger -> info("LinkedCellContainer destructed");
    }
};

TEST_F(LinkedCellContainerTest, Analytical1) {
    test_logger->info("LinkedCellContainer - Two body analytical test");
    container2D->addParticle(Particle(std::array<double, 3>{7, 7.5, 0.5}, std::array<double, 3>{-1, 0, 0}, 1));
    container2D->addParticle(Particle(std::array<double, 3>{8, 7.5, 0.5}, std::array<double, 3>{1, 0, 0}, 1));
    container2D->simulate(1, 1e-6, "", "", 10, false, true);
    double expectedVel = sqrt(80 * (pow(1.0 / 3, 6) - pow(1.0 / 3, 12)) + 4) / 2;
    EXPECT_NEAR(-expectedVel, container2D->getParticles()[0].getV()[0], 1e-4);
    EXPECT_NEAR(0, container2D->getParticles()[0].getV()[1], 1e-12);
    EXPECT_NEAR(0, container2D->getParticles()[0].getV()[2], 1e-12);
    EXPECT_NEAR(expectedVel, container2D->getParticles()[1].getV()[0], 1e-4);
    EXPECT_NEAR(0, container2D->getParticles()[1].getV()[1], 1e-12);
    EXPECT_NEAR(0, container2D->getParticles()[1].getV()[2], 1e-12);
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Two body analytical test failed");
    } else {
        test_logger->info("LinkedCellContainer - Two body analytical test passed");
    }

}
