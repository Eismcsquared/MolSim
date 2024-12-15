#include <gtest/gtest.h>
#include <string>
#include "inputReader/XMLReader.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"
#include "container/DirectSumContainer.h"
#include "container/LinkedCellContainer.h"
#include "body/Sphere.h"

class XMLReaderTest: public ::testing::Test {
protected:
    std::vector<Particle> particles;
    std::string inputFile;
    std::unique_ptr<Simulation> simulation;
    const double pi = 3.141592653589793;
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

    void SetUp() override {
        particles = std::vector<Particle>();
    }
};

// Test whether the input file for the assignment 1 is correctly parsed.
TEST_F(XMLReaderTest, Assigment1Input) {
    test_logger->info("XMLReader - Assignment 1 input test");
    inputFile = "../tests/test_cases/assignment1.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    EXPECT_EQ(4, simulation->getContainer()->getParticleNumber());
    std::vector<Particle> ref_p;
    ref_p.emplace_back(std::array<double, 3>{0, 0, 0}, std::array<double, 3>{0, 0, 0}, 1);
    ref_p.emplace_back(std::array<double, 3>{0, 1, 0}, std::array<double, 3>{-1, 0, 0}, 3e-6);
    ref_p.emplace_back(std::array<double, 3>{0, 5.36, 0}, std::array<double, 3>{-0.425, 0, 0}, 9.55e-4);
    ref_p.emplace_back(std::array<double, 3>{34.75, 0, 0}, std::array<double, 3>{0, 0.0296, 0}, 1e-14);

    std::unique_ptr<Force> f = std::make_unique<GravitationalForce>();

    DirectSumContainer ref(ref_p, f);
    EXPECT_EQ(ref, *(simulation->getContainer()));
    EXPECT_NEAR(1000, simulation->getEndTime(), 1e-12);
    EXPECT_NEAR(0.014, simulation->getDeltaT(), 1e-12);
    EXPECT_EQ("MD_vtk", simulation->getOutputFile());
    EXPECT_EQ("vtu", simulation->getOutputFormat());
    EXPECT_EQ(10, simulation->getOutputFrequency());
    EXPECT_FLOAT_EQ(0, simulation->getContainer()->getG());
    EXPECT_FALSE(simulation->getContainer()->getThermostat());
    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Assignment 1 input test failed\n\n");
    } else {
        test_logger->info("XMLReader - Assignment 1 input test passed\n\n");
    }
}

// Test whether the input file for the assignment 2 is correctly parsed. The Brownian motion is set to 0 for test purpose.
TEST_F(XMLReaderTest, Assigment2Input) {
    test_logger->info("XMLReader - Assignment 2 input test");
    inputFile = "../tests/test_cases/assignment2.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    EXPECT_EQ(384, simulation->getContainer()->getParticleNumber());
    std::vector<Particle> ref_p;
    Cuboid c1(std::array<double, 3>{0, 0, 0.5}, std::array<double, 3>{0, 0, 0}, std::array<unsigned int, 3>{40, 8, 1}, 1, pow(2, 1.0 / 6), 0, 2);
    Cuboid c2(std::array<double, 3>{15, 15, 0.5}, std::array<double, 3>{0, -10, 0}, std::array<unsigned int, 3>{8, 8, 1}, 1, pow(2, 1.0 / 6), 0, 2);
    c1.createParticles(ref_p);
    c2.createParticles(ref_p);
    std::unique_ptr<Force> f = std::make_unique<LennardJonesForce>();
    DirectSumContainer ref(ref_p, f);
    EXPECT_EQ(ref, *(simulation->getContainer()));
    EXPECT_NEAR(5, simulation->getEndTime(), 1e-12);
    EXPECT_NEAR(0.0002, simulation->getDeltaT(), 1e-12);
    EXPECT_EQ("a2", simulation->getOutputFile());
    EXPECT_EQ("xyz", simulation->getOutputFormat());
    EXPECT_EQ(15, simulation->getOutputFrequency());
    EXPECT_FLOAT_EQ(0, simulation->getContainer()->getG());
    EXPECT_FALSE(simulation->getContainer()->getThermostat());
    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Assignment 2 input test failed\n\n");
    } else {
        test_logger->info("XMLReader - Assignment 2 input test passed\n\n");
    }

}

// Test whether the input file for the assignment 3 is correctly parsed. The Brownian motion is set to 0 for test purpose.
TEST_F(XMLReaderTest, Assigment3Input) {
    test_logger->info("XMLReader - Assignment 3 input test");
    inputFile = "../tests/test_cases/assignment3.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    EXPECT_EQ(2400, simulation->getContainer()->getParticleNumber());
    std::vector<Particle> ref_p;
    Cuboid c1(std::array<double, 3>{20, 20, 0.5}, std::array<double, 3>{0, 0, 0}, std::array<unsigned int, 3>{100, 20, 1}, 1, pow(2, 1.0 / 6), 0, 2);
    Cuboid c2(std::array<double, 3>{70, 60, 0.5}, std::array<double, 3>{0, -10, 0}, std::array<unsigned int, 3>{20, 20, 1}, 1, pow(2, 1.0 / 6), 0, 2);
    c1.createParticles(ref_p);
    c2.createParticles(ref_p);
    std::unique_ptr<Force> f = std::make_unique<LennardJonesForce>();
    LinkedCellContainer ref(ref_p, f, std::array<double, 3>{180, 90, 1}, 3, std::array<BoundaryCondition, 6>{OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW});
    EXPECT_EQ(ref, dynamic_cast<LinkedCellContainer&>(*(simulation->getContainer())));
    EXPECT_NEAR(20, simulation->getEndTime(), 1e-12);
    EXPECT_NEAR(0.0005, simulation->getDeltaT(), 1e-12);
    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Assignment 3 input test failed\n\n");
    } else {
        test_logger->info("XMLReader - Assignment 3 input test passed\n\n");
    }
}

// Test whether the input file for the fallig drop simulation is correctly parsed. The Brownian motion is set to 0 for test purpose.
TEST_F(XMLReaderTest, FallingDropInput) {
    test_logger->info("XMLReader - Falling drop input test");
    inputFile = "../tests/test_cases/falling_drop.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    std::vector<Particle> ref_p;
    Sphere s(std::array<double, 3>{60, 25, 0.5}, std::array<double, 3>{0, -10, 0}, 15, 1, pow(2, 1.0 / 6), 0, 2);
    s.createParticles(ref_p);
    std::unique_ptr<Force> f = std::make_unique<LennardJonesForce>();
    LinkedCellContainer ref(ref_p, f, std::array<double, 3>{120, 50, 1}, 3, std::array<BoundaryCondition, 6>{REFLECTING, REFLECTING, REFLECTING, REFLECTING, OUTFLOW, OUTFLOW});
    EXPECT_EQ(ref, dynamic_cast<LinkedCellContainer&>(*(simulation->getContainer())));
    EXPECT_NEAR(10, simulation->getEndTime(), 1e-12);
    EXPECT_NEAR(0.00005, simulation->getDeltaT(), 1e-12);
    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Falling drop input test failed\n\n");
    } else {
        test_logger->info("XMLReader - Falling drop input test passed\n\n");
    }
}

// Test whether Brownian motion is initialized correctly.
TEST_F(XMLReaderTest, BrownianMotion) {
    test_logger->info("XMLReader - Brownian motion reading test");
    inputFile = "../tests/test_cases/large.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    long double mean = 0;
    long double meanSquire = 0;
    EXPECT_EQ(10000, particles.size());
    for(Particle &p: particles) {
        mean += ArrayUtils::L2Norm(p.getV());
        meanSquire += pow(ArrayUtils::L2Norm(p.getV()), 2);
    }
    mean /= particles.size();
    meanSquire /= particles.size();
    EXPECT_NEAR(0.1 * sqrt(pi / 2), mean, 2e-3) << "Wrong mean velocity of Brownian motion. Expected: "
    << 0.1 * sqrt(pi / 2) << ", but got " << mean;
    EXPECT_NEAR(0.02, meanSquire, 5e-4) << "Wrong mean squared velocity of Brownian motion. Expected: "
    << 0.02 << ", but got " << meanSquire;
    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Brownian motion reading test failed\n\n");
    } else {
        test_logger->info("XMLReader - Brownian motion reading test passed\n\n");
    }
}

// Test whether the xml reader read the thermostat and gravitational acceleration correctly

TEST_F(XMLReaderTest, Assignment4Input) {
    test_logger->info("XMLReader - Assignment 4 input test");
    inputFile = "../tests/test_cases/assignment4.xml";
    simulation = XMLReader::readXML(particles, inputFile);
    EXPECT_TRUE(simulation->getContainer()->getThermostat());
    EXPECT_EQ(1000, simulation->getContainer()->getThermostat()->getPeriode());
    EXPECT_FLOAT_EQ(40, simulation->getContainer()->getThermostat()->getTargetT());
    EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), simulation->getContainer()->getThermostat()->getMaxDelta());
    EXPECT_EQ(2, simulation->getContainer()->getThermostat()->getDimension());
    EXPECT_FLOAT_EQ(-12.44, simulation->getContainer()->getG());

    // Test whether the parameters epsilon, sigma and the type of particles are read correcly as well.
    for (Particle &p: particles) {
        EXPECT_FLOAT_EQ(1, p.getEpsilon());
        if (std::abs(p.getM() - 1) < 1e-12) {
            EXPECT_FLOAT_EQ(1, p.getSigma());
            EXPECT_EQ(1, p.getType());
        } else {
            EXPECT_FLOAT_EQ(0.9412, p.getSigma());
            EXPECT_EQ(2, p.getType());
        }
    }

    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Assignment 4 input test failed\n\n");
    } else {
        test_logger->info("XMLReader - Assignment 4 input test passed\n\n");
    }
}

// Test the generation of Brownian motion by temperature

TEST_F(XMLReaderTest, InitialTemperature) {
    test_logger->info("XMLReader - Initial temperature test");
    inputFile = "../tests/test_cases/large_with_thermostat.xml";
    particles.reserve(20000);
    simulation = XMLReader::readXML(particles, inputFile);

    EXPECT_TRUE(simulation->getContainer()->getThermostat());
    EXPECT_EQ(1000, simulation->getContainer()->getThermostat()->getPeriode());
    EXPECT_FLOAT_EQ(60, simulation->getContainer()->getThermostat()->getTargetT());
    EXPECT_FLOAT_EQ(0.95, simulation->getContainer()->getThermostat()->getMaxDelta());
    EXPECT_EQ(3, simulation->getContainer()->getThermostat()->getDimension());

    // two different masses

    double meanSquare1 = 0;
    double meanSquare2 = 0;

    for (Particle &p: simulation->getContainer()->getParticles()) {
        if (std::abs(p.getM() - 1) < 1e-12) {
            meanSquare1 += pow(ArrayUtils::L2Norm(p.getV()), 2);
        } else {
            meanSquare2 += pow(ArrayUtils::L2Norm(p.getV()), 2);
        }
    }

    meanSquare1 /= 1e5;
    meanSquare2 /= 1e5;

    EXPECT_NEAR(120, meanSquare1, 1);
    EXPECT_NEAR(60, meanSquare2, .5);

    if(::testing::Test::HasFailure()) {
        test_logger->info("XMLReader - Initial temperature test failed\n\n");
    } else {
        test_logger->info("XMLReader - Initial temperature test passed\n\n");
    }
}



