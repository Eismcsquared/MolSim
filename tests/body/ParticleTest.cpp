#include <spdlog/spdlog.h>
#include <gtest/gtest.h>
#include "body/Particle.h"

class ParticleTest: public ::testing::Test {
protected:
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");
};

// Test whether particles are created with the correct identity.
TEST_F(ParticleTest, ID) {
    test_logger->info("Particle - ID test");
    Particle p0(1);
    Particle p1({0, 0, 0}, {0, 0, 0}, 1);
    Particle p2({1, 1, 1}, {1, 1, 1}, 2, 0, 5, 2);
    int id0 = p0.getId();
    EXPECT_EQ(id0 + 1, p1.getId());
    EXPECT_EQ(id0 + 2, p2.getId());
    Particle p1_copy = p1;
    EXPECT_EQ(id0 + 1, p1_copy.getId());
    EXPECT_TRUE(p1.is(p1_copy));
    Particle p3(2);
    EXPECT_EQ(id0 + 3, p3.getId());

    if (::testing::Test::HasFailure()) {
        test_logger->info("Particle - ID test failed\n\n");
    } else {
        test_logger->info("Particle - ID test passed\n\n");
    }
}

// Test that stationary particles do not move.
TEST_F(ParticleTest, Stationary) {
    test_logger->info("Particle - Stationary particle test");
    Particle p({1, -3, 9}, {4, 1, -1}, 1, 0, 5, 1, true);
    p.setF({10, -6, 3});
    p.updateX(1);
    p.updateV(1);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(p.getX() - std::array<double, 3>{1, -3, 9}), 1e-12);
    EXPECT_NEAR(0, ArrayUtils::L2Norm(p.getV() - std::array<double, 3>{4, 1, -1}), 1e-12);

    if (::testing::Test::HasFailure()) {
        test_logger->info("Particle - Stationary particle test failed\n\n");
    } else {
        test_logger->info("Particle - Stationary particle test passed\n\n");
    }
}
