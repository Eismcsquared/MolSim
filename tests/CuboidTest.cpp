#include <vector>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "Logger.h"
#include "FileReader.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "LennardJonesForce.h"

class CuboidTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;
    Cuboid cuboid = Cuboid({0, 10, 0}, {0, 0, 0}, {2, 4, 3}, 1, 1, 0);
    char* testfile = const_cast<char*>("../tests/test_cases/two_cuboid.txt");


    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        pc = new ParticleContainer(particles, new LennardJonesForce());
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("Particle Container created");
    }

    void TearDown() override {
        test_logger->info("Particle Container deleted\n\n");
        delete pc;
    }
};

// Test whether the createParticle method of Cuboid create the correct particles that are contained in the cuboid.
TEST_F(CuboidTest, CreateParticle) {
    test_logger->info("Cuboid - create particles test");
    std::vector<Particle> ps;
    cuboid.createParticles(ps);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 3; ++k) {
                Particle ref = Particle({static_cast<double>(i), static_cast<double>(10 + j), static_cast<double>(k)}, {0, 0, 0}, 1);
                ASSERT_EQ(ref, ps[12 * i + 3 * j + k]) << "Cuboid - create particles test failed";
            }
        }
    }
    test_logger->info("Cuboid - create particles test passed");
}

// Test whether the extended input format is recognized correctly.
TEST_F(CuboidTest, ReadCuboids) {
    test_logger->info("Cuboid - read file test");
    std::vector<Particle> ps;
    std::array<double, 3> x = {0, 0, 0};
    std::array<double, 3> v = {0, 0, 0};
    ps.emplace_back(x, v, 1);
    v = {1.5, -1.5, 0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 1; ++k) {
                x = {static_cast<double>(6 + i), static_cast<double>(j), static_cast<double>(k)};
                ps.emplace_back(x, v, 1);
            }
        }
    }
    ParticleContainer reference(ps, new LennardJonesForce());
    ASSERT_EQ(reference, *pc) << "Cuboid - read file test failed";
    test_logger->info("Cuboid - read file test passed");
}

// Test whether the method addCuboid correctly adds a cuboid to the container
TEST_F(CuboidTest, AddCuboid) {
    test_logger->info("Cuboid - add cuboid test");
    std::vector<Particle> ps;
    pc->addCuboid(cuboid);
    std::array<double, 3> x = {0, 0, 0};
    std::array<double, 3> v = {0, 0, 0};
    ps.emplace_back(x, v, 1);
    v = {1.5, -1.5, 0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 1; ++k) {
                x = {static_cast<double>(6 + i), static_cast<double>(j), static_cast<double>(k)};
                ps.emplace_back(x, v, 1);
            }
        }
    }
    v = {0, 0, 0};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 3; ++k) {
                x = {static_cast<double>(i), static_cast<double>(10 + j), static_cast<double>(k)};
                ps.emplace_back(x, v, 1);
            }
        }
    }
    ParticleContainer reference(ps, new LennardJonesForce());
    ASSERT_EQ(reference, *pc) << "Cuboid - read file test failed";
    test_logger->info("Cuboid - read file test filed");
}

// Test whether the method createParticle of Cuboid generate particles with Brownian motion with the specified average speed
TEST_F(CuboidTest, BrownianMotion) {
    test_logger->info("Cuboid - Brownian motion test");
    Cuboid c = Cuboid({0, 0, 0}, {0, 0, 0}, {100, 100, 100}, 1, 1, 1); // 1000000 particles should produce statistically significant results
    std::vector<Particle> ps;
    c.createParticles(ps);
    double averageSpeedSquare = 0;
    for (const Particle& p : ps) {
        averageSpeedSquare += pow(ArrayUtils::L2Norm(p.getV()), 2);
    }
    averageSpeedSquare /= 1000000;
    std::cout << averageSpeedSquare;
    ASSERT_FALSE(std::abs(averageSpeedSquare - 3) > 0.01) << "Cuboid - Brownian motion test failed";
    test_logger->info("Cuboid - Brownian motion test passed");
}