#include <gtest/gtest.h>
#include "Particle.h"
#include "ParticleContainer.h"
#include "GravitationalForce.h"
#include <vector>
#include <iostream>
#include "FileReader.h"
#include "utils/ArrayUtils.h"
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/fmt/ostr.h>
#include <fstream>
#include <sstream>
#include "Logger.h"

class SimpleTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;
    char* testfile = const_cast<char*>("../tests/test_cases/two_body.txt");


    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        pc = new ParticleContainer(particles, new GravitationalForce());
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("Particle Container created");
    }

    void TearDown() override {
        test_logger->info("delete Particle Container\n\n");
        delete pc;
    }
};

//Test that the FileReader still works after modifying it to read cuboids
TEST_F(SimpleTest, ReadFile) {
    test_logger->info("Read file test");
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 2; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    ParticleContainer reference(ref_vec, new GravitationalForce());
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


// Test whether the method addParticle of ParticleContainer works when adding a single particle
TEST_F(SimpleTest, AddParticle1) {
    test_logger->info("Add particle test1");
    Particle p = Particle({2,0,0}, {0,0,0}, 1, 0);
    pc->addParticle(p);
    std::vector<Particle> ref_vec;
    for (int i = 0; i < 3; ++i) {
        std::array<double, 3> x = {static_cast<double>(i), 0, 0};
        std::array<double, 3> v = {0, 0, 0};
        ref_vec.emplace_back(x, v, 1);
    }
    ParticleContainer reference(ref_vec, new GravitationalForce());
    if(*pc == reference) {
        test_logger->info("Add particle test1 passed");
    } else {
        test_logger->error("Add particle test1 failed");
        test_logger->error("Expected: " + reference.toString());
        test_logger->error("But got: " + pc->toString());
    }
    ASSERT_EQ(*pc, reference);
}
