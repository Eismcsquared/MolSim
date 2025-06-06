#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <fstream>
#include <algorithm>
#include "statistics/Statistics.h"
#include "container/LinkedCellContainer.h"
#include "force/LennardJonesForce.h"

class StatisticsTest: public ::testing::Test {
protected:
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");
    std::vector<Particle> particles;
    std::unique_ptr<LinkedCellContainer> container;

    void SetUp() override {
        particles.emplace_back(std::array<double, 3>{1, 3, 6}, std::array<double, 3>{3, -2, 6}, 1);
        particles.emplace_back(std::array<double, 3>{32, 12, 44}, std::array<double, 3>{4, 3, 2}, 1);
        // not in domain.
        particles.emplace_back(std::array<double, 3>{55, 10, 28}, std::array<double, 3>{12, 2, 5}, 1);
        particles.emplace_back(std::array<double, 3>{8, 42, 9}, std::array<double, 3>{-1, -4, 3}, 1);
        particles.emplace_back(std::array<double, 3>{21, 43, 49}, std::array<double, 3>{7, 9, 6}, 1);
        particles.emplace_back(std::array<double, 3>{48, 12, 13}, std::array<double, 3>{-2, 5, -7}, 1);
        particles.emplace_back(std::array<double, 3>{49, 31, 37}, std::array<double, 3>{3, 1, 0}, 1);
        particles.emplace_back(std::array<double, 3>{2, 23, 19}, std::array<double, 3>{-10, -4, 6}, 1);
        particles.emplace_back(std::array<double, 3>{2, 37, 8}, std::array<double, 3>{3, -5, -1}, 1);
        particles.emplace_back(std::array<double, 3>{6, 33, 16}, std::array<double, 3>{-10, -4, -8}, 1);
        // stationary particle should not be considered.
        particles.emplace_back(std::array<double, 3>{16, 13, 22}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
        particles.emplace_back(std::array<double, 3>{42, 28, 11}, std::array<double, 3>{-2, -7, 2}, 1);
        particles.emplace_back(std::array<double, 3>{32, 21, 22}, std::array<double, 3>{-8, 1, -10}, 1);
        particles.emplace_back(std::array<double, 3>{36, 33, 42}, std::array<double, 3>{10, -3, 4}, 1);
        particles.emplace_back(std::array<double, 3>{44, 3, 22}, std::array<double, 3>{5, 8, 8}, 1);
        particles.emplace_back(std::array<double, 3>{28, 34, 8}, std::array<double, 3>{-6, 0, -3}, 1);
        particles.emplace_back(std::array<double, 3>{33, 48, 6}, std::array<double, 3>{4, -6, 9}, 1);
        std::unique_ptr<Force> force = std::make_unique<LennardJonesForce>();
        container = std::make_unique<LinkedCellContainer>(particles, force, std::array<double, 3>{50, 50, 50}, 3,
                                                          std::array<BoundaryCondition, 6>{REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING});
    }
};

// Test the computation of the density and the velocity profile.
TEST_F(StatisticsTest, Computation) {
    test_logger->info("Statistics - Computation test");
    Statistics s1("", 10, 0, 50, 10);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(s1.density(container->getParticles()) - std::vector<double>{
            0.6, 0.4, 0, 0, 0.2, 0.2, 0.6, 0.2, 0.4, 0.4
    }), 1e-12);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(s1.velocity(container->getParticles()) - std::vector<double>{
            -11.0 / 3, -4, 0, 0, 9, 0, -2.0 / 3, -3, 0.5, 3
    }), 1e-12);
    Statistics s2("", 10, 0, 50, 10, Y, Z);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(s2.density(container->getParticles()) - std::vector<double>{
            0.4, 0, 0.4, 0, 0.4, 0.2, 0.8, 0.2, 0.4, 0.2
    }), 1e-12);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(s2.velocity(container->getParticles()) - std::vector<double>{
            7, 0, -2.5, 0, -2, 2, -1.75, -1, 4.5, 9
    }), 1e-12);

    if (::testing::Test::HasFailure()) {
        test_logger->info("Statistics - Computation test failed\n\n");
    } else {
        test_logger->info("Statistics - Computation test passed\n\n");
    }
}

// Test whether the statistics is correctly written into the specified output file
TEST_F(StatisticsTest, Output) {
    test_logger->info("Statistics - Output test");
    Statistics s("s", 10, 0, 50, 10);
    s.saveStatistics(container->getParticles(), 30);
    std::array<double, 10> density{0.6, 0.4, 0, 0, 0.2, 0.2, 0.6, 0.2, 0.4, 0.4};
    std::array<double, 10> velocity{-11.0 / 3, -4, 0, 0, 9, 0, -2.0 / 3, -3, 0.5, 3};
    std::string fileName = "s_0030.csv";
    std::ifstream inFile(fileName);
    for (int i = 0; i < 10; ++i) {
        std::string line;
        std::getline(inFile, line);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream dataStream(line);
        double x, rho, v;
        dataStream >> x;
        dataStream >> rho;
        dataStream >> v;
        EXPECT_NEAR(5 * i + 2.5, x, 1e-12);
        EXPECT_NEAR(density[i], rho, 1e-12);
        EXPECT_NEAR(velocity[i], v, 1e-12);
    }
    std::remove(fileName.c_str());

    if (::testing::Test::HasFailure()) {
        test_logger->info("Statistics - Output test failed\n\n");
    } else {
        test_logger->info("Statistics - Output test passed\n\n");
    }
}

// Test whether the statistics unit is correctly integrated into the simulation.
TEST_F(StatisticsTest, Integration) {
    test_logger->info("Statistics - Integration test");

    std::shared_ptr<Statistics> s = std::make_shared<Statistics>("s", 10, 0, 50, 10);
    container->simulate(0.01, 1e-3, "", "", 100, true, 0, s);

    std::string fileName = "s_0010.csv";
    std::ifstream inFile(fileName);
    for (int i = 0; i < 10; ++i) {
        std::string line;
        std::getline(inFile, line);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream dataStream(line);
        double x;
        dataStream >> x;
        EXPECT_NEAR(5 * i + 2.5, x, 1e-12);
    }
    std::remove(fileName.c_str());

    if (::testing::Test::HasFailure()) {
        test_logger->info("Statistics - Integration test failed\n\n");
    } else {
        test_logger->info("Statistics - Integration test passed\n\n");
    }
}
