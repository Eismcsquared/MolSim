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

char* testfile = "../input/simpletest.txt";

class PartilceTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;

    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        GravitationalForce f;
        pc = new ParticleContainer(particles, 0, 1000, 0.014, f, ".vtu");
    }

    void TearDown() override {
        std::cout << "delete Particle Container " << "\n";
        delete pc;
    }
};

TEST_F(PartilceTest, Addparticle) {
    Particle p = Particle({1,0,0}, {0,0,0}, 1, 0);
    int size_before = pc->getParticleSize();
    pc->addParticle(p);
    int size_after = pc->getParticleSize();
    ASSERT_EQ(size_before + 1, size_after) ;
}

TEST_F(PartilceTest, Addparticle_2) {
    std::cout << pc->getParticleSize() << "\n";
    int size_before = pc->getParticleSize();
    for(int i = 0 ; i < 10; i++) {
        Particle p = Particle({(double)i +2,0,0}, {0,0,0}, 1, 0);
        pc->addParticle(p);
    }
    int size_after = pc->getParticleSize();
    ASSERT_EQ(size_before + 10, size_after) ;
}