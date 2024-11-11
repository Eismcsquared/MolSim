#include <gtest/gtest.h>
#include "Particle.h"
#include "ParticleContainer.h"
#include "GravitationalForce.h"
#include <vector>
#include <iostream>
#include "FileReader.h"
#include "utils/ArrayUtils.h"
#include <string>
//#include "outputWriter/XYZWriter.h"
//#include "outputWriter/VTKWriter.h"

char* testfile = "test.txt";

class PartilceTest : public ::testing::Test {
protected:

    std::vector<Particle> particles;
    ParticleContainer *particle_container;

    void SetUp() override {
        FileReader fileReader;
        ASSERT_NO_THROW(fileReader.readFile(particles, testfile)) << "Error reading file";
        particle_container = new ParticleContainer(particles, 0, 1000, 0.014, new GravitationalForce(), ".vtu");
    }

    void TearDown() override {
        std::cout << "delete Particle Container " << "\n";
        //delete particle_container;
    }
};

TEST_F(PartilceTest, Test1) {
    std::cout << "Test1 : " <<  "\n";
    ASSERT_EQ(particles.size(),3 ) << "Particle size not same ";


}

