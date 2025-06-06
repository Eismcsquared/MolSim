#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/fmt/ostr.h>
#include "force/LennardJonesForce.h"
#include "force/Force.h"
#include "body/Particle.h"
#include "container/DirectSumContainer.h"
#include "force/GravitationalForce.h"

class LennardJonesForceTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    std::unique_ptr<ParticleContainer> pc;
    std::unique_ptr<Force> f;
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");


    void SetUp() override {
        f = std::make_unique<LennardJonesForce>();
        pc = std::make_unique<DirectSumContainer>(particles, f);
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("ParticleContainer created");
    }

    void TearDown() override {
        test_logger->info("ParticleContainer deleted\n\n");
    }
};

// Test the Calculation of the Lennard Jones force for a two body system
TEST_F(LennardJonesForceTest, ForceCalculation1) {
    test_logger->info("Lennard-Jones force - two-body test");
    pc->addParticle(Particle({0, 0, 0}, {0, 0, 0}, 1));
    pc->addParticle(Particle({1, 0, 0}, {0, 0, 0}, 1));
    pc->updateF();
    Particle ref_p1({0, 0, 0}, {0, 0, 0}, 1);
    ref_p1.setF({-120, 0, 0}); // expected force of particle 1
    Particle ref_p2({1, 0, 0}, {0, 0, 0}, 1);
    ref_p2.setF({120, 0, 0}); // expected force of particle 2
    if (pc->getParticles()[0] == ref_p1 && pc->getParticles()[1] == ref_p2) {
        test_logger->info("Lennard-Jones force - two-body test passed");
    } else {
        test_logger->error("Lennard-Jones force - two-body test failed");
    }
    ASSERT_EQ(ref_p1, pc->getParticles()[0]);
    ASSERT_EQ(ref_p2, pc->getParticles()[1]);
}

// Test the Calculation of the Lennard Jones force for multiple particles
TEST_F(LennardJonesForceTest, ForceCalculation2) {
    test_logger->info("Lennard-Jones force - four-body test");
    pc->addParticle(Particle({0, 0, 0}, {0, 0, 0}, 1));
    pc->addParticle(Particle({1, 0, 0}, {0, 0, 0}, 1));
    pc->addParticle(Particle({0, 1, 0}, {0, 0, 0}, 1));
    pc->addParticle(Particle({1, 1, 0}, {0, 0, 0}, 1));
    pc->updateF();
    Particle ref_p1({0, 0, 0}, {0, 0, 0}, 1);
    ref_p1.setF({-915.0 / 8, -915.0 / 8, 0}); // expected force of particle 1
    Particle ref_p2({1, 0, 0}, {0, 0, 0}, 1);
    ref_p2.setF({915.0 / 8, -915.0 / 8, 0}); // expected force of particle 2
    Particle ref_p3({0, 1, 0}, {0, 0, 0}, 1);
    ref_p3.setF({-915.0 / 8, 915.0 / 8, 0}); // expected force of particle 3
    Particle ref_p4({1, 1, 0}, {0, 0, 0}, 1);
    ref_p4.setF({915.0 / 8, 915.0 / 8, 0}); // expected force of particle 4
    if (pc->getParticles()[0] == ref_p1 && pc->getParticles()[1] == ref_p2 &&
        pc->getParticles()[2] == ref_p3 && pc->getParticles()[3] == ref_p4) {
        test_logger->info("Lennard-Jones force - four-body test passed");
    } else {
        test_logger->error("Lennard-Jones force - four-body test failed");
    }
    ASSERT_EQ(ref_p1, pc->getParticles()[0]);
    ASSERT_EQ(ref_p2, pc->getParticles()[1]);
    ASSERT_EQ(ref_p3, pc->getParticles()[2]);
    ASSERT_EQ(ref_p4, pc->getParticles()[3]);
}

// Test whether the Lorentz-Berthelot mixing rule is correctly applied in the force calculation.
TEST_F(LennardJonesForceTest, LorentzBerthelotMixingRule) {
    test_logger->info("Lennard-Jones force - Lorentz-Berthelot mixing rule test");
    pc->addParticle(Particle({0, 0, 0}, {0, 0, 0}, 1, 0, 5, 1));
    pc->addParticle(Particle({1, 0, 0}, {0, 0, 0}, 1, 0, 3, 2));
    pc->updateF();
    double expected = 24 * sqrt(15) * (2 * pow(1.5, 12) - pow(1.5, 6));
    EXPECT_LE(ArrayUtils::L2Norm(pc->getParticles()[0].getF() - std::array<double, 3>{-expected, 0, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(pc->getParticles()[1].getF() - std::array<double, 3>{expected, 0, 0}), 1e-12);
    if (::testing::Test::HasFailure()) {
        test_logger->info("Lennard-Jones force - Lorentz-Berthelot mixing rule test failed");
    } else {
        test_logger->info("Lennard-Jones force - Lorentz-Berthelot mixing rule test passed");
    }
}
