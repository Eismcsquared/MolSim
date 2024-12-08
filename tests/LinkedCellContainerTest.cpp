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

//Test whether cells are instantiated with correct size and number, especially if the cell size does not divide the domain size.
TEST_F(LinkedCellContainerTest, CellInstantiation) {
    test_logger->info("LinkedCellContainer - Cell instantiation test");
    EXPECT_EQ(5, container2D->getNCells()[0]);
    EXPECT_EQ(5, container2D->getNCells()[1]);
    EXPECT_EQ(1, container2D->getNCells()[2]);
    for (auto c: container2D->getCells()) {
        EXPECT_EQ(3, c.getSize()[0]);
        EXPECT_EQ(3, c.getSize()[1]);
        EXPECT_EQ(1, c.getSize()[2]);
    }
    EXPECT_EQ(147, container2D->getCells().size());

    EXPECT_EQ(33, container3D->getNCells()[0]);
    EXPECT_EQ(33, container3D->getNCells()[1]);
    EXPECT_EQ(16, container3D->getNCells()[2]);
    for (auto c: container3D->getCells()) {
        EXPECT_FLOAT_EQ(100.0 / 33, c.getSize()[0]);
        EXPECT_FLOAT_EQ(100.0 / 33, c.getSize()[1]);
        EXPECT_FLOAT_EQ(50.0 / 16, c.getSize()[2]);
    }
    EXPECT_EQ(22050, container3D->getCells().size());
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Cell instantiation test failed");
    } else {
        test_logger->info("LinkedCellContainer - Cell instantiation test passed");
    }
}

// Test the transformation between 1D and 3D indices.
TEST_F(LinkedCellContainerTest, Indices) {
    test_logger->info("LinkedCellContainer - Get indices test");
    EXPECT_EQ(0, container2D->get1DIndex({-1, -1, -1}));
    EXPECT_EQ(13583, container3D->get1DIndex({2, 2, 10}));
    EXPECT_EQ(64, container2D->getCellIndex({1, 3.5, 0.4}));
    EXPECT_EQ(-1, container3D->getCellIndex({105, 50, 40}));
    std::array<int, 3> index3D = container2D->get3DIndex(28);
    EXPECT_EQ(-1, index3D[0]);
    EXPECT_EQ(3, index3D[1]);
    EXPECT_EQ(-1, index3D[2]);
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Get indices test failed");
    } else {
        test_logger->info("LinkedCellContainer - Get indices test passed");
    }
}

TEST_F(LinkedCellContainerTest, Neighbours) {
    test_logger->info("LinkedCellContainer - Neighbour indices test");
    // inner cell
    std::set<int> ref1 = {11636, 11638, 11671, 11672, 11673, 11601, 11602, 11603,
                          12862, 12861, 12863, 12896, 12897, 12898, 12826, 12827, 12828,
                          10412, 10411, 10413, 10446, 10447, 10448, 10376, 10377, 10378};
    EXPECT_EQ(ref1, container3D->getCells()[11637].getNeighbours());
    // cell at boundary
    std::set<int> ref2 = {1523, 1525, 1558, 1559, 1560, 1488, 1489, 1490,
                          2749, 2748, 2750, 2783, 2784, 2785, 2713, 2714, 2715};
    EXPECT_EQ(ref2, container3D->getCells()[1524].getNeighbours());
    // cell at vertex
    std::set<int> ref3 = {7387, 7421, 7422, 8611, 8612, 8646, 8647, 6161, 6162, 6196, 6197};
    EXPECT_EQ(ref3, container3D->getCells()[7386].getNeighbours());
    // cell at corner
    std::set<int> ref4 = {20787, 20753, 20752, 19563, 19562, 19528, 19527};
    EXPECT_EQ(ref4, container3D->getCells()[20788].getNeighbours());
    // halo cells
    std::set<int> empty{};
    EXPECT_EQ(empty, container3D->getCells()[0].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[337].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[1237].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[12244].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[5740].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[5424].getNeighbours());
}

// Test whether the method removeFrom works correctly
TEST_F(LinkedCellContainerTest, RemoveHalo) {
    test_logger->info("LinkedCellContainer - Remove halo particles test");
    // In halo cells of all 6 direction plus two positions in domain.
    std::array<std::array<double, 3>, 8> pos = {
            std::array<double, 3>{-1, 40, 46},
            std::array<double, 3>{102, 75, 30},
            std::array<double, 3>{60, -2, 40},
            std::array<double, 3>{46, 101, 28},
            std::array<double, 3>{70, 70, -1},
            std::array<double, 3>{63, 49, 50.5},
            std::array<double, 3>{63, 49, 38},
            std::array<double, 3>{74, 23, 19}};
    for (int i = 0; i < 8; ++i) {
        container3D->addParticle(Particle(pos[i], {0, 0, 0}, 1));
    }
    EXPECT_EQ(8, container3D->getParticleNumber());
    for (int i = 0; i < 8; ++i) {
        EXPECT_EQ(1, container3D->getCells()[container3D->getCellIndex(pos[i])].getParticleIndices().size());
        EXPECT_TRUE(container3D->getParticles()[i].isInDomain());
    }

    for (int i = 0; i < 6; ++i) {
        // remove particles in halo cells direction by direction.
        container3D->removeFromHalo(static_cast<Direction>(i));
        for (int j = 0; j <= i; ++j) {
            EXPECT_EQ(0, container3D->getCells()[container3D->getCellIndex(pos[j])].getParticleIndices().size());
            EXPECT_FALSE(container3D->getParticles()[j].isInDomain());
        }
        for (int j = i + 1; j < 8; ++j) {
            EXPECT_EQ(1, container3D->getCells()[container3D->getCellIndex(pos[j])].getParticleIndices().size());
            EXPECT_TRUE(container3D->getParticles()[j].isInDomain());
        }
        // Also check the getParticleNumber method.
        EXPECT_EQ(7 - i, container3D->getParticleNumber());
    }
}

// Test whether particles more than cutoff apart indeed does not affect each other
TEST_F(LinkedCellContainerTest, ForceCalculation) {
    test_logger->info("LinkedCellContainer - Force calculation test");
    std::array<std::array<double, 3>, 8> pos = {
            std::array<double, 3>{10, 10, 10},
            std::array<double, 3>{11, 11, 11},
            std::array<double, 3>{13, 9, 9} // "for away" form the first two
    };

    // updateF called in addParticle, test the calculation
    container3D->addParticle(Particle(pos[0], {0, 0, 0}, 1));
    container3D->addParticle(Particle(pos[1], {0, 0, 0}, 1));
    double expectedForce = 120 * (pow(1 / sqrt(3), 7) - 2 * pow(1 / sqrt(3), 13));
    for (int i = 0; i < 3; ++i) {
        EXPECT_FLOAT_EQ(expectedForce / sqrt(3), container3D->getParticles()[0].getF()[i]);
        EXPECT_FLOAT_EQ(-expectedForce / sqrt(3), container3D->getParticles()[1].getF()[i]);
    }
    container3D->addParticle(Particle(pos[2], {0, 0, 0}, 1));
    // expected no interaction between the third particle and the previous two, since they are too far away.
    for (int i = 0; i < 3; ++i) {
        EXPECT_FLOAT_EQ(expectedForce / sqrt(3), container3D->getParticles()[0].getF()[i]);
        EXPECT_FLOAT_EQ(-expectedForce / sqrt(3), container3D->getParticles()[1].getF()[i]);
        EXPECT_FLOAT_EQ(0, container3D->getParticles()[2].getF()[i]);
    }
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Force calculation test failed");
    } else {
        test_logger->info("LinkedCellContainer - Force calculation test passed");
    }
}

// Test the outflow boundary condition
TEST_F(LinkedCellContainerTest, Outflow) {
    test_logger->info("LinkedCellContainer - Outflow test");
    container2D->addParticle(Particle({1, 1, 0.5}, {-1, 0, 0}, 1)); // outflow at t=1
    container2D->addParticle(Particle({13, 13, 0.5}, {0, 1, 0}, 1)); // outflow at t=2
    container2D->addParticle(Particle({13, 1, 0.5}, {sqrt(2), sqrt(2), 0}, 1)); // outflow at t=sqrt(2)
    EXPECT_EQ(3, container2D->getParticleNumber());
    container2D->simulate(0.99, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(3, container2D->getParticleNumber());
    container2D->simulate(0.02, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(2, container2D->getParticleNumber());
    container2D->simulate(0.40, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(2, container2D->getParticleNumber());
    container2D->simulate(0.01, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(1, container2D->getParticleNumber());
    container2D->simulate(0.57, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(1, container2D->getParticleNumber());
    container2D->simulate(0.02, 1e-3, "", "", 10, false, true);
    EXPECT_EQ(0, container2D->getParticleNumber());
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Outflow test failed");
    } else {
        test_logger->info("LinkedCellContainer - Outflow test passed");
    }
}

TEST_F(LinkedCellContainerTest, Reflecting) {
    test_logger->info("LinkedCellContainer - Reflecting test");
    container3D->addParticle(Particle({1, 1, 10}, {-1, 0, 0}, 1)); //reflected at t=1
    container3D->addParticle(Particle({98, 98, 10}, {1, 1, 0}, 1)); //reflected at t=2
    container3D->addParticle(Particle({50, 1, 48}, {2, -1, 1}, 1)); //reflect at t=1, outflow at t=2
    container3D->simulate(1.5, 1e-2, "", "", 10, false, true);
    EXPECT_NEAR(0.5, container3D->getParticles()[0].getX()[0], 2e-2);
    EXPECT_FLOAT_EQ(1, container3D->getParticles()[0].getV()[0]);
    EXPECT_NEAR(53, container3D->getParticles()[2].getX()[0], 2e-2);
    EXPECT_NEAR(0.5, container3D->getParticles()[2].getX()[1], 2e-2);
    EXPECT_NEAR(49.5, container3D->getParticles()[2].getX()[2], 2e-2);
    EXPECT_FLOAT_EQ(1, container3D->getParticles()[2].getV()[1]);
    EXPECT_EQ(3, container3D->getParticleNumber());
    container3D->simulate(1, 1e-2, "", "", 10, false, true);
    EXPECT_NEAR(99.5, container3D->getParticles()[1].getX()[0], 2e-2);
    EXPECT_NEAR(99.5, container3D->getParticles()[1].getX()[1], 2e-2);
    EXPECT_FLOAT_EQ(-1, container3D->getParticles()[1].getV()[0]);
    EXPECT_FLOAT_EQ(-1, container3D->getParticles()[1].getV()[1]);
    EXPECT_EQ(2, container3D->getParticleNumber());
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Reflecting test failed");
    } else {
        test_logger->info("LinkedCellContainer - Reflecting test passed");
    }
}

// Test whether the potential is correctly cut of at the cutoff radius on a two body system. The idea is, after both particles
// leaving the potential of each other, their kinetic energy should be the sum of their initial kinetic energy plus the potential
// they pass throw.
TEST_F(LinkedCellContainerTest, Analytical) {
    test_logger->info("LinkedCellContainer - Two body analytical test");
    container2D->addParticle(Particle(std::array<double, 3>{7, 7.5, 0.5}, std::array<double, 3>{0, 0, 0}, 1, 0.25, 1));
    container2D->addParticle(Particle(std::array<double, 3>{8, 7.5, 0.5}, std::array<double, 3>{0, 0, 0}, 1, 0.25, 1));
    container2D->simulate(25, 1e-3, "", "", 10, false, true);
    double expectedVel = sqrt((pow(1.0 / 3, 6) - pow(1.0 / 3, 12)));
    EXPECT_TRUE(ArrayUtils::L2Norm(container2D->getParticles()[0].getX() - container2D->getParticles()[1].getX()) > 3);
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

