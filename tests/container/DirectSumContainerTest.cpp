#include <vector>
#include <string>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "inputReader/FileReader.h"
#include "body/Particle.h"
#include "container/DirectSumContainer.h"
#include "force/GravitationalForce.h"

class DirectSumContainerTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    std::unique_ptr<Force> f ;
    std::unique_ptr<DirectSumContainer> pc;
    std::string testfile = "../tests/test_cases/two_body.txt";
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");


    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        f = std::make_unique<GravitationalForce>();
        pc = std::make_unique<DirectSumContainer>(particles, f);
        pc->updateF();
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("ParticleContainer created");
    }

    void TearDown() override {
        test_logger->info("ParticleContainer deleted\n\n");
    }
};

//Test that the FileReader still works after modifying it to read cuboids. especially if additional white spaces are appended
// to the actual data
TEST_F(DirectSumContainerTest, ReadFile) {
    test_logger->info("DirectSumContainer - Read file test");
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 2; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    std::unique_ptr<Force> ref_f = std::make_unique<GravitationalForce>();
    DirectSumContainer reference(ref_vec, ref_f);
    reference.updateF();
    if(*pc == reference) {
        test_logger->info("DirectSumContainer - Read file test passed");
    } else {
        test_logger->error("DirectSumContainer - Read file test failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Test whether the calculation of gravitational forces for a simple example works (which is done in the constructor of DirectSumContainer)
TEST_F(DirectSumContainerTest, Gravitation) {
    test_logger->info("DirectSumContainer - Gravitation test");
    Particle ref_p1({0, 0, 0}, {0, 0, 0}, 1);
    ref_p1.setF({1, 0, 0}); // expected force of particle 1
    Particle ref_p2({1, 0, 0}, {0, 0, 0}, 1);
    ref_p2.setF({-1, 0, 0}); // expected force of particle 2
    if (pc->getParticles()[0] == ref_p1 && pc->getParticles()[1] == ref_p2) {
        test_logger->info("DirectSumContainer - Gravitation test passed");
    } else {
        test_logger->error("DirectSumContainer - Gravitation test failed");
    }
    ASSERT_EQ(ref_p1, pc->getParticles()[0]);
    ASSERT_EQ(ref_p2, pc->getParticles()[1]);
}

// Test whether the method addParticle of DirectSumContainer works when adding a single particle
TEST_F(DirectSumContainerTest, AddParticle1) {
    test_logger->info("DirectSumContainer - Add particle test 1");
    Particle p = Particle({2,0,0}, {0,0,0}, 1, 0);
    pc->addParticle(p);
    pc->updateF();
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 3; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }

    std::unique_ptr<Force> ref_f = std::make_unique<GravitationalForce>();
    DirectSumContainer reference(ref_vec, ref_f);
    reference.updateF();
    if(*pc == reference) {
        test_logger->info("DirectSumContainer - Add particle test 1 passed");
    } else {
        test_logger->error("DirectSumContainer - Add particle test 1 failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Test whether the method addParticle of DirectSumContainer works when adding multiple particles
TEST_F(DirectSumContainerTest, AddParticle2) {
    test_logger->info("DirectSumContainer - Add particle test 2");
    for(int i = 1 ; i <= 10; i++) {
        Particle p = Particle({(double) i + 1, 0, 0}, {0, 0, 0}, 1, 0);
        pc->addParticle(p);
    }
    pc->updateF();
    std::vector<Particle> ref_vec;
    for (int i = 0; i <= 11; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }

    std::unique_ptr<Force> ref_f = std::make_unique<GravitationalForce>();
    DirectSumContainer reference(ref_vec, ref_f);
    reference.updateF();
    if(*pc == reference) {
        test_logger->info("DirectSumContainer - Add particle test 2 passed");
    } else {
        test_logger->error("DirectSumContainer - Add particle test 2 failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Compare the simulation of this simple example with its analytic solution, one can show that for x = x_2 - x_1 it hols:
// arccos(sqrt(x)) + sqrt(x) * sqrt(1 - x) = 2 * t => x = 1 / 2 (hence x_1 = 0.25, x_2 = 0.75 due to symmetry) if t = pi / 8 + 1 / 4
TEST_F(DirectSumContainerTest, Analytical) {
    test_logger->info("DirectSumContainer - Analytical solution test");
    double pi = 3.14159265358979323846;
    double end_t = pi / 8 + 0.25;
    double delta_t = 1e-6;
    pc->simulate(end_t, delta_t, "", "vtu", 10, false);
    std::vector<Particle> ref_vec;
    std::array<double, 3> x_1 = {0.25, 0, 0}; // expected position of particle 1
    std::array<double, 3> x_2 = {0.75, 0, 0}; // expected position of particle 2
    std::array<double, 3> v_1 = {1, 0, 0}; // expected velocity of particle 1
    std::array<double, 3> v_2 = {-1, 0, 0}; // expected velocity of particle 2
    ref_vec.emplace_back(x_1, v_1, 1);
    ref_vec.emplace_back(x_2, v_2, 1);

    std::unique_ptr<Force> ref_f = std::make_unique<GravitationalForce>();
    DirectSumContainer reference(ref_vec, ref_f);
    reference.updateF();
    if (!(*pc == reference)) {
        test_logger->error("DirectSumContainer - Analytical solution test failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
    test_logger->info("DirectSumContainer - Analytical solution test passed");
}

// Analytical test: The motion of a planet in the gravitational field generated by a stationary sun. We set
// M = m = 1, G = 1, r_0 = (1, 0, 0), v_0 = (0, 1/2, 0)
// The orbit is given as
// r(\phi) = 1 / ( 4 - 3 * cos(\phi))
// from which the identity follows
// |r| + |r - (6/7, 0, 0)| = 8/7
TEST_F(DirectSumContainerTest, AnalyticalOneBody) {
    test_logger->info("DirectSumContainer - One body analytical test");
    std::vector<Particle> p{};
    std::unique_ptr<Force> gravitation = std::make_unique<GravitationalForce>();
    DirectSumContainer container(p, gravitation);
    // The sun
    container.addParticle(Particle({0, 0, 0}, {0, 0, 0}, 1, 0, 5, 1, true));
    // The planet
    container.addParticle(Particle({1, 0, 0}, {0, 0.5, 0}, 1));

    for (int i = 0; i < 10; ++i) {
        container.simulate(i, 1e-5, "", "", 10, false);
        EXPECT_NEAR(8.0 / 7, ArrayUtils::L2Norm(container.getParticles()[1].getX()) + ArrayUtils::L2Norm(container.getParticles()[1].getX() - std::array<double, 3>{6.0 / 7, 0, 0}), 1e-5);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(container.getParticles()[0].getX()), 1e-12);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(container.getParticles()[0].getV()), 1e-12);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(container.getParticles()[0].getF()), 1e-12);
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("DirectSumContainer - One body analytical test failed");
    } else {
        test_logger->info("DirectSumContainer - One body analytical test passed");
    }
}