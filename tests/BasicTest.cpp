#include <vector>
#include <string>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "Logger.h"
#include "FileReader.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "GravitationalForce.h"

class BasicTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;
    char* testfile = const_cast<char*>("../tests/test_cases/two_body.txt");


    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        pc = new ParticleContainer(particles, new GravitationalForce(),true);
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("Particle Container created");
    }

    void TearDown() override {
        test_logger->info("Particle Container deleted\n\n");
        delete pc;
    }
};

//Test that the FileReader still works after modifying it to read cuboids. especially if additional white spaces are appended
// to the actual data
TEST_F(BasicTest, ReadFile) {
    test_logger->info("Read file test");
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 2; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    ParticleContainer reference(ref_vec, new GravitationalForce(),true);
    if(*pc == reference) {
        test_logger->info("Read file test passed");
    } else {
        test_logger->error("Read file test failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Test whether the calculation of gravitational forces for a simple example works (which is done in the constructor of ParticleContainer)
TEST_F(BasicTest, Gravitation) {
    test_logger->info("Gravitation test");
    Particle ref_p1({0, 0, 0}, {0, 0, 0}, 1);
    ref_p1.setF({1, 0, 0}); // expected force of particle 1
    Particle ref_p2({1, 0, 0}, {0, 0, 0}, 1);
    ref_p2.setF({-1, 0, 0}); // expected force of particle 2
    if (pc->getParticles()[0] == ref_p1 && pc->getParticles()[1] == ref_p2) {
        test_logger->info("Gravitation test passed");
    } else {
        test_logger->error("Gravitation test failed");
    }
    ASSERT_EQ(ref_p1, pc->getParticles()[0]);
    ASSERT_EQ(ref_p2, pc->getParticles()[1]);
}

// Test whether the method addParticle of ParticleContainer works when adding a single particle
TEST_F(BasicTest, AddParticle1) {
    test_logger->info("Add particle test 1");
    Particle p = Particle({2,0,0}, {0,0,0}, 1, 0);
    pc->addParticle(p);
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 3; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    ParticleContainer reference(ref_vec, new GravitationalForce(),true);
    if(*pc == reference) {
        test_logger->info("Add particle test 1 passed");
    } else {
        test_logger->error("Add particle test 1 failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Test whether the method addParticle of ParticleContainer works when adding multiple particles
TEST_F(BasicTest, AddParticle2) {
    test_logger->info("Add particle test 2");
    for(int i = 1 ; i <= 10; i++) {
        Particle p = Particle({(double) i + 1, 0, 0}, {0, 0, 0}, 1, 0);
        pc->addParticle(p);
    }
    std::vector<Particle> ref_vec;
    for (int i = 0; i <= 11; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    ParticleContainer reference(ref_vec, new GravitationalForce(),true);
    if(*pc == reference) {
        test_logger->info("Add particle test 2 passed");
    } else {
        test_logger->error("Add particle test 2 failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}

// Compare the simulation of this simple example with its analytic solution, one can show that for x = x_2 - x_1 it hols:
// arccos(sqrt(x)) + sqrt(x) * sqrt(1 - x) = 2 * t => x = 1 / 2 (hence x_1 = 0.25, x_2 = 0.75 due to symmetry) if t = pi / 8 + 1 / 4
TEST_F(BasicTest, Analytic) {
    test_logger->info("Analytic solution test");
    double pi = 3.14159265358979323846;
    double end_t = pi / 8 + 0.25;
    double delta_t = 1e-6;
    for (int i = 0; i < (int) (end_t / delta_t); ++i) {
        pc->updateX(delta_t);
        pc->updateF();
        pc->updateV(delta_t);
    }
    std::vector<Particle> ref_vec;
    std::array<double, 3> x_1 = {0.25, 0, 0}; // expected position of particle 1
    std::array<double, 3> x_2 = {0.75, 0, 0}; // expected position of particle 2
    std::array<double, 3> v_1 = {1, 0, 0}; // expected velocity of particle 1
    std::array<double, 3> v_2 = {-1, 0, 0}; // expected velocity of particle 2
    ref_vec.emplace_back(x_1, v_1, 1);
    ref_vec.emplace_back(x_2, v_2, 1);
    ParticleContainer reference(ref_vec, new GravitationalForce(),true);
    if (!(*pc == reference)) {
        test_logger->error("Analytic solution test failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
    test_logger->info("Analytic solution test passed");
}

