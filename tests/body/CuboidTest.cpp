#include <vector>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "inputReader/FileReader.h"
#include "utils/ArrayUtils.h"
#include "body/Particle.h"
#include "container/DirectSumContainer.h"
#include "force/LennardJonesForce.h"

class CuboidTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    std::unique_ptr<DirectSumContainer> pc;
    std::unique_ptr<Force> f;
    Cuboid cuboid = Cuboid({0, 10, 0}, {0, 0, 0}, {2, 4, 3}, 1, 1, 0, 3);
    std::string testfile = "../tests/test_cases/two_cuboid.txt";
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        f = std::make_unique<LennardJonesForce>();
        pc = std::make_unique<DirectSumContainer>(particles, f);
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("ParticleContainer created");
    }

    void TearDown() override {
        test_logger->info("ParticleContainer deleted\n\n");
    }
};

// Test whether the createParticle method of Cuboid create the correct particles that are contained in the cuboid.
TEST_F(CuboidTest, CreateParticle) {
    test_logger->info("Cuboid - Create particles test");
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
    test_logger->info("Cuboid - Create particles test passed");
}

// Test whether the extended input format is recognized correctly.
TEST_F(CuboidTest, ReadCuboids) {
    test_logger->info("Cuboid - Read file test");
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

    std::unique_ptr<Force> ref_f = std::make_unique<LennardJonesForce>();
    DirectSumContainer reference(ps, ref_f);
    ASSERT_EQ(reference, *pc) << "Cuboid - Read file test failed";
    test_logger->info("Cuboid - Read file test passed");
}

// Test whether the method addCuboid correctly adds a cuboid to the container
TEST_F(CuboidTest, AddCuboid) {
    test_logger->info("Cuboid - Add cuboid test");
    std::vector<Particle> ps;
    pc->addCluster(cuboid);
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
    std::unique_ptr<Force> ref_f = std::make_unique<LennardJonesForce>();
    DirectSumContainer reference(ps, ref_f);
    ASSERT_EQ(reference, *pc) << "Cuboid - Read file test failed";
    test_logger->info("Cuboid - Read file test filed");
}

// Test whether the method createParticle of Cuboid generate particles with Brownian motion with the specified average speed
// Note that since the average speed is not easy to compute, we compute the average of squared speed instead.
// This can be computed as the integral over (v_x^2 + v_y^2 + v_z^2) * exp(-(v_x^2 + v_y^2 + v_z^2) / 2) / sqrt(2 * pi)^3,
// assuming passing 1 as the last parameter to Cuboid. This integral is just a sum of Gaussian integrals and can be evaluated
// analytically. The result is 3.
TEST_F(CuboidTest, BrownianMotion) {
    test_logger->info("Cuboid - Brownian motion test");
    Cuboid c = Cuboid({0, 0, 0}, {0, 0, 0}, {100, 100, 100}, 1, 1, 1, 3); // 1000000 particles should produce statistically significant results
    std::vector<Particle> ps;
    c.createParticles(ps);
    double averageSpeedSquare = 0;
    for (const Particle& p : ps) {
        averageSpeedSquare += pow(ArrayUtils::L2Norm(p.getV()), 2);
    }
    averageSpeedSquare /= 1000000;
    ASSERT_FALSE(std::abs(averageSpeedSquare - 3) > 0.01) << "Cuboid - Brownian motion test failed";
    test_logger->info("Cuboid - Brownian motion test passed");
}

// Test the creation of stationary cuboids.
TEST_F(CuboidTest, Stationary) {
    test_logger->info("Cuboid - Stationary cuboid test");
    Cuboid c({0, 0, 0}, {0, 0, 0}, {3, 2, 1}, 1, 1, 0, 3, 0, 5, 1, true);
    std::vector<Particle> ps;
    c.createParticles(ps);
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ps[i].getX() - std::array<double, 3>{static_cast<double>(i >> 1), static_cast<double>(i % 2), 0}), 1e-12);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ps[i].getV() - std::array<double, 3>{0, 0, 0}), 1e-12);
        EXPECT_TRUE(ps[i].isStationary());
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("Cuboid - Stationary cuboid test failed");
    } else {
        test_logger->info("Cuboid - Stationary cuboid test passed");
    }
}