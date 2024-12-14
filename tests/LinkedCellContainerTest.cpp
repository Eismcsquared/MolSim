#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include "container/LinkedCellContainer.h"
#include "force/LennardJonesForce.h"
#include "Logger.h"

class LinkedCellContainerTest: public ::testing::Test {
protected:
    std::vector<Particle> p1;
    std::vector<Particle> p2;
    std::vector<Particle> p3;
    std::vector<Particle> p4;
    std::unique_ptr<LinkedCellContainer> container2D;
    std::unique_ptr<LinkedCellContainer> container3D;
    std::unique_ptr<LinkedCellContainer> container3DPeriodicX;
    std::unique_ptr<LinkedCellContainer> container3DPeriodicAll;

    void SetUp() override {

        std::unique_ptr<Force> f1 = std::make_unique<LennardJonesForce>();
        std::unique_ptr<Force> f2 = std::make_unique<LennardJonesForce>();
        std::unique_ptr<Force> f3 = std::make_unique<LennardJonesForce>();
        std::unique_ptr<Force> f4 = std::make_unique<LennardJonesForce>();

        container2D = std::make_unique<LinkedCellContainer>(p1, f1, std::array<double, 3>{15, 15, 1}, 3,
                                          std::array<BoundaryCondition, 6>{OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW});
        container3D = std::make_unique<LinkedCellContainer>(p2, f2, std::array<double, 3>{100, 100, 50}, 3,
                                          std::array<BoundaryCondition, 6>{REFLECTING, REFLECTING, REFLECTING, REFLECTING, OUTFLOW, OUTFLOW});
        container3DPeriodicX = std::make_unique<LinkedCellContainer>(p3, f3, std::array<double, 3>{100, 100, 50}, 3,
                                                                     std::array<BoundaryCondition, 6>{PERIODIC, PERIODIC, REFLECTING, REFLECTING, OUTFLOW, OUTFLOW});
        container3DPeriodicAll = std::make_unique<LinkedCellContainer>(p4, f4, std::array<double, 3>{100, 100, 50}, 3,
                                                                       std::array<BoundaryCondition, 6>{PERIODIC, PERIODIC, PERIODIC, PERIODIC, PERIODIC, PERIODIC});
        test_logger -> info("LinkedCellContainer created");
    }

    void TearDown() override {
        test_logger -> info("LinkedCellContainer destructed");
    }
};

//Test whether cells are instantiated with correct size and number, especially if the cell size does not divide the domain size.
TEST_F(LinkedCellContainerTest, CellInstantiation) {
    test_logger->info("LinkedCellContainer - Cell instantiation test");
    EXPECT_TRUE((std::array<int, 3>{5, 5, 1}) == container2D->getNCells());
    for (const auto& c: container2D->getCells()) {
        EXPECT_LE(ArrayUtils::L2Norm(c.getSize() - std::array<double, 3>{3, 3, 1}), 1e-12);
    }
    EXPECT_EQ(147, container2D->getCells().size());

    EXPECT_TRUE((std::array<int, 3>{33, 33, 16}) == container3D->getNCells());
    for (const auto& c: container3D->getCells()) {
        EXPECT_LE(ArrayUtils::L2Norm(c.getSize() - std::array<double, 3>{100.0 / 33, 100.0 / 33, 50.0 / 16}), 1e-12);
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
    EXPECT_TRUE((std::array<int, 3>{-1, 3, -1}) == index3D);
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Get indices test failed");
    } else {
        test_logger->info("LinkedCellContainer - Get indices test passed");
    }
}

// Test whether neighbouring cells are computed correctly.
TEST_F(LinkedCellContainerTest, Neighbours) {
    test_logger->info("LinkedCellContainer - Neighbour indices test");
    // inner cell
    std::set<int> ref1 = {11636, 11638, 11671, 11672, 11673, 11601, 11602, 11603,
                          12862, 12861, 12863, 12896, 12897, 12898, 12826, 12827, 12828,
                          10412, 10411, 10413, 10446, 10447, 10448, 10376, 10377, 10378};
    EXPECT_EQ(ref1, container3D->getCells()[11637].getNeighbours());
    EXPECT_EQ(ref1, container3DPeriodicX->getCells()[11637].getNeighbours());
    // cell at boundary
    std::set<int> ref2 = {1523, 1525, 1558, 1559, 1560, 1488, 1489, 1490,
                          2749, 2748, 2750, 2783, 2784, 2785, 2713, 2714, 2715};
    EXPECT_EQ(ref2, container3D->getCells()[1524].getNeighbours());
    EXPECT_EQ(ref2, container3DPeriodicX->getCells()[1524].getNeighbours());
    // cell at vertex
    std::set<int> ref3 = {7387, 7421, 7422, 8611, 8612, 8646, 8647, 6161, 6162, 6196, 6197};
    std::set<int> ref3Periodic = {7387, 7418, 7421, 7422, 7453, 8611, 8612, 8643, 8646, 8647, 8678, 6161, 6162, 6193, 6196, 6197, 6228};
    EXPECT_EQ(ref3, container3D->getCells()[7386].getNeighbours());
    EXPECT_EQ(ref3Periodic, container3DPeriodicX->getCells()[7386].getNeighbours());
    // cell at corner
    std::set<int> ref4 = {20787, 20753, 20752, 19563, 19562, 19528, 19527};
    std::set<int> ref4Periodic = {20787, 20756, 20753, 20752, 20721, 19563, 19562, 19531, 19528, 19527, 19496};
    EXPECT_EQ(ref4, container3D->getCells()[20788].getNeighbours());
    EXPECT_EQ(ref4Periodic, container3DPeriodicX->getCells()[20788].getNeighbours());
    // halo cells
    std::set<int> empty{};
    EXPECT_EQ(empty, container3D->getCells()[0].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[337].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[1237].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[12244].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[5740].getNeighbours());
    EXPECT_EQ(empty, container3D->getCells()[5424].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[0].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[337].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[1237].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[12244].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[5740].getNeighbours());
    EXPECT_EQ(empty, container3DPeriodicX->getCells()[5424].getNeighbours());
    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Neighbour indices test failed");
    } else {
        test_logger->info("LinkedCellContainer - Neighbour indices test passed");
    }
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
            std::array<double, 3>{13, 9, 9} // "far away" form the first two
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

// Test the correctness of the force calculation in case of periodic boundary condition.
TEST_F(LinkedCellContainerTest, ForceCalculationPeriodic) {
    test_logger->info("LinkedCellContainer - Periodic force calculation test");
    std::array<std::array<double, 3>, 8> pos = {
            std::array<double, 3>{1, 2, 2},
            std::array<double, 3>{99, 2, 2},
            std::array<double, 3>{99, 2, 4},
            std::array<double, 3>{99, 4, 2},
            std::array<double, 3>{1, 99, 2}, // no interaction with the first particle since the y-boundaries are not periodic
            std::array<double, 3>{1, 2, 49}, // no interaction with the first particle since the z-boundaries are not periodic
    };

    auto force = [](double r) {
        return 120 * (pow(r, -7) - 2 * pow(r, -13));
    };

    container3DPeriodicX->addParticle(Particle(pos[0], {0, 0, 0}, 1));
    container3DPeriodicX->addParticle(Particle(pos[1], {0, 0, 0}, 1));

    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[0].getF() -
                 std::array<double, 3>{
                    -force(2),
                    0,
                    0
                 }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[1].getF() -
                std::array<double, 3>{
                    force(2),
                    0,
                    0
                }),
              1e-12);

    container3DPeriodicX->addParticle(Particle(pos[2], {0, 0, 0}, 1));

    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[0].getF() -
                std::array<double, 3>{
                    -force(2) - force(sqrt(8)) / sqrt(2),
                    0,
                    force(sqrt(8)) / sqrt(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[1].getF() -
                std::array<double, 3>{
                    force(2),
                    0,
                    force(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[2].getF() -
                std::array<double, 3>{
                    force(sqrt(8)) / sqrt(2),
                    0,
                    -force(2) - force(sqrt(8)) / sqrt(2)
                }),
              1e-12);

    container3DPeriodicX->addParticle(Particle(pos[3], {0, 0, 0}, 1));


    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[0].getF() -
                std::array<double, 3>{
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[1].getF() -
                std::array<double, 3>{
                    force(2),
                    force(2),
                    force(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[2].getF() -
                std::array<double, 3>{
                    force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2),
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[3].getF() -
                std::array<double, 3>{
                    force(sqrt(8)) / sqrt(2),
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2)
                }),
              1e-12);


    container3DPeriodicX->addParticle(Particle(pos[4], {0, 0, 0}, 1));
    container3DPeriodicX->addParticle(Particle(pos[5], {0, 0, 0}, 1));

    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[0].getF() -
                std::array<double, 3>{
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[1].getF() -
                std::array<double, 3>{
                    force(2),
                    force(2),
                    force(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[2].getF() -
                std::array<double, 3>{
                    force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2),
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2)
                }),
              1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[3].getF() -
                std::array<double, 3>{
                    force(sqrt(8)) / sqrt(2),
                    -force(2) - 2 * force(sqrt(8)) / sqrt(2),
                    force(sqrt(8)) / sqrt(2)
                }),
              1e-12);

    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[4].getF() -
            std::array<double, 3>{
                    0,
                    0,
                    0
                }),
              1e-12);

    EXPECT_LE(ArrayUtils::L2Norm(container3DPeriodicX->getParticles()[4].getF() -
                std::array<double, 3>{
                    0,
                    0,
                    0
                }),
              1e-12);

    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Periodic force calculation test failed");
    } else {
        test_logger->info("LinkedCellContainer - Periodic force calculation test passed");
    }
}

// Test the outflow boundary condition
TEST_F(LinkedCellContainerTest, Outflow) {
    test_logger->info("LinkedCellContainer - Outflow boundary test");
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
        test_logger->info("LinkedCellContainer - Outflow boundary test failed");
    } else {
        test_logger->info("LinkedCellContainer - Outflow boundary test passed");
    }
}

// Test the reflecting boundary condition
TEST_F(LinkedCellContainerTest, Reflecting) {
    test_logger->info("LinkedCellContainer - Reflecting boundary test");
    container3D->addParticle(Particle({1, 1, 10}, {-1, 0, 0}, 1)); //reflected at t=1
    container3D->addParticle(Particle({98, 98, 10}, {1, 1, 0}, 1)); //reflected at t=2
    container3D->addParticle(Particle({50, 1, 48}, {2, -1, 1}, 1)); //reflect at t=1, outflow at t=2
    container3D->simulate(1.5, 1e-1, "", "", 10, false, true);

    EXPECT_EQ(3, container3D->getParticleNumber());
    EXPECT_EQ(Particle({0.5, 1, 10}, {1, 0, 0}, 1), container3D->getParticles()[0]);
    EXPECT_EQ(Particle({99.5, 99.5, 10}, {1, 1, 0}, 1), container3D->getParticles()[1]);
    EXPECT_EQ(Particle({53, 0.5, 49.5}, {2, 1, 1}, 1), container3D->getParticles()[2]);

    container3D->simulate(1, 1e-1, "", "", 10, false, true);

    EXPECT_EQ(2, container3D->getParticleNumber());
    EXPECT_FALSE(container3D->getParticles()[2].isInDomain());
    EXPECT_EQ(Particle({1.5, 1, 10}, {1, 0, 0}, 1), container3D->getParticles()[0]);
    EXPECT_EQ(Particle({99.5, 99.5, 10}, {-1, -1, 0}, 1), container3D->getParticles()[1]);

    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Reflecting boundary test failed");
    } else {
        test_logger->info("LinkedCellContainer - Reflecting boundary test passed");
    }
}

// Test the periodic boundary condition
TEST_F(LinkedCellContainerTest, Periodic) {
    test_logger->info("LinkedCellContainer - Periodic boundary test");
    container3DPeriodicX->addParticle(Particle({1, 1, 10}, {-1, 0, 0}, 1)); // periodic at t=1
    container3DPeriodicX->addParticle(Particle({98, 98, 10}, {1, 1, 0}, 1)); // periodic + reflecting at t=2
    container3DPeriodicX->addParticle(Particle({1, 1, 48}, {-1, 0, 2}, 1)); // periodic + outflow at t=1

    container3DPeriodicX->simulate(1.5, 1e-1, "", "", 10, false, true);
    EXPECT_EQ(2, container3DPeriodicX->getParticleNumber());
    EXPECT_FALSE(container3DPeriodicX->getParticles()[2].isInDomain());
    EXPECT_EQ(Particle({99.5, 1, 10}, {-1, 0, 0}, 1), container3DPeriodicX->getParticles()[0]);
    EXPECT_EQ(Particle({99.5, 99.5, 10}, {1, 1, 0}, 1), container3DPeriodicX->getParticles()[1]);
    container3DPeriodicX->simulate(1, 1e-1, "", "", 10, false, true);
    EXPECT_EQ(2, container3DPeriodicX->getParticleNumber());
    EXPECT_EQ(Particle({98.5, 1, 10}, {-1, 0, 0}, 1), container3DPeriodicX->getParticles()[0]);
    EXPECT_EQ(Particle({0.5, 99.5, 10}, {1, -1, 0}, 1), container3DPeriodicX->getParticles()[1]);

    // What happens at "periodic corners"?

    container3DPeriodicAll->addParticle(Particle({1, 1, 1}, {-1, -1, -1}, 1));
    container3DPeriodicAll->simulate(2, 1e-1, "", "", 10, false, true);
    EXPECT_EQ(Particle({99, 99, 49}, {-1, -1, -1}, 1), container3DPeriodicAll->getParticles()[0]);

    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Periodic boundary test failed");
    } else {
        test_logger->info("LinkedCellContainer - Periodic boundary test passed");
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

    EXPECT_LE(ArrayUtils::L2Norm(container2D->getParticles()[0].getV() -
                std::array<double, 3>{
                    -expectedVel,
                    0,
                    0
                }),
              1e-4);

    EXPECT_LE(ArrayUtils::L2Norm(container2D->getParticles()[1].getV() -
                 std::array<double, 3>{
                     expectedVel,
                     0,
                     0
                 }),
               1e-4);

    if (::testing::Test::HasFailure()) {
        test_logger->info("LinkedCellContainer - Two body analytical test failed");
    } else {
        test_logger->info("LinkedCellContainer - Two body analytical test passed");
    }
}

