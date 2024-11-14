#include <gtest/gtest.h>
#include "Particle.h"
#include "ParticleContainer.h"
#include "GravitationalForce.h"
#include <vector>
#include <iostream>
#include "FileReader.h"
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "Logger.h"



// Test whether the simulations in assignment 1 still work as before, reference data originates from previous simulations
class Assignment1Test : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;
    char* testfile = const_cast<char*>("../tests/test_cases/assignment1.txt");
    FileReader fileReader;

    void SetUp() override {
        fileReader.readFile(particles, testfile);
        pc = new ParticleContainer(particles, new GravitationalForce(), true);
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("Particle Container created");
    }

    void TearDown() override {
        test_logger->info("Particle Container deleted\n\n");
        delete pc;
    }
};

// Test the calculation of the gravitational forces
TEST_F(Assignment1Test, Assignment1_updateF) {
    test_logger->info("Assignment 1 - update force test");
    Particle ref_p1({0, 0, 0}, {0, 0, 0}, 1);
    ref_p1.setF({8.28114e-18, 3.6241e-05, 0}); // expected force of the sun
    Particle ref_p2({0, 1, 0}, {-1, 0, 0}, 3e-06);
    ref_p2.setF({2.48126e-23, -2.99985e-06, 0}); // expected force of the earth
    Particle ref_p3({0, 5.36, 0}, {-0.425, 0, 0}, 0.000955);
    ref_p3.setF({7.63443e-21, -3.32411e-05, 0}); // expected force of the jupiter
    Particle ref_p4({34.75, 0, 0}, {0, 0.0296, 0}, 1e-14);
    ref_p4.setF({-8.2888e-18, 1.17828e-21, 0}); // expected force of the Halley comet
    if (pc->getParticles()[0] == ref_p1 && pc->getParticles()[1] == ref_p2 &&
    pc->getParticles()[2] == ref_p3 && pc->getParticles()[3] == ref_p4) {
        test_logger->info("Assignment 1 - update force test passed");
    } else {
        test_logger->error("Assignment 1 - update force test failed");
    }
    ASSERT_EQ(ref_p1, pc->getParticles()[0]);
    ASSERT_EQ(ref_p2, pc->getParticles()[1]);
    ASSERT_EQ(ref_p3, pc->getParticles()[2]);
    ASSERT_EQ(ref_p4, pc->getParticles()[3]);
}

// Test the state of the particles after several iterations
TEST_F(Assignment1Test, Simulation_simple) {
    test_logger->info("Assignment 1 - simple simulation test");
    for(int iteration = 0 ; iteration <= 3; ++iteration) {
        pc->updateX(0.014);
        pc->updateF();
        pc->updateV(0.014);
    }
    std::vector<Particle> ref_p;
    char *ref_file = const_cast<char*>("../tests/Answer_Ref/Ans_simulation_simple.txt");
    fileReader.readFile(ref_p, ref_file);
    ParticleContainer reference(ref_p, new GravitationalForce(), true);

    if(!(reference == *pc)) {
        test_logger->error("Assignment 1 - simple simulation test failed");
        std::cout << reference.toString() << std::endl;
        std::cout << pc->toString() << std::endl;
    }
    ASSERT_EQ(reference, *pc) << "Assignment 1 - simple simulation test failed";
    test_logger->info("Assignment 1 - simple simulation test passed");
}

// Test the state of the particles after the entire simulation
TEST_F(Assignment1Test, Complex_simulation) {
    test_logger->info("Assignment 1 - complex simulation test");
    int max_iteration = (int) ((1000 - 0) / 0.014);
    for(int iteration = 0 ; iteration <= max_iteration; ++iteration) {
        pc->updateX(0.014);
        pc->updateF();
        pc->updateV(0.014);
    }
    std::vector<Particle> ref_p;
    char *ref_file = const_cast<char*>("../tests/Answer_Ref/Ans_simulation_complex.txt");
    fileReader.readFile(ref_p, ref_file);
    ParticleContainer reference(ref_p, new GravitationalForce(),true);

    if(!(reference == *pc)) {
        test_logger->error("Assignment 1 - complex simulation test failed");
        std::cout << reference.toString() << std::endl;
        std::cout << pc->toString() << std::endl;
    }
    ASSERT_EQ(reference, *pc) << "Assignment 1 - complex simulation test failed";
    test_logger->info("Assignment 1 - complex simulation test passed");
}