#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include "container/LinkedCellContainer.h"
#include "force/LennardJonesForce.h"
#include "inputReader/StateReader.h"
#include "inputReader/XMLReader.h"
#include "Simulation.h"

class ParallelizationTest: public ::testing::Test {
protected:
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

};

// Run a simulation and test whether the parallelized versions produce the same result as the sequential one.
TEST_F(ParallelizationTest, Simulation) {
    test_logger->info("Parallelization - Simulation test");

    std::string input1 = "../tests/test_cases/parallelization1.xml";
    std::string input2 = "../tests/test_cases/parallelization2.xml";
    std::string inputReference = "../tests/test_cases/reference.txt";

    std::vector<Particle> ref{};
    StateReader::loadState(ref, inputReference);
    EXPECT_EQ(648, ref.size());

    // parallelization strategy 0
    std::vector<Particle> p1{};
    std::unique_ptr<Simulation> simulation1 = XMLReader::readXML(p1, input1);
    // parallelization strategy 1
    std::vector<Particle> p2{};
    std::unique_ptr<Simulation> simulation2 = XMLReader::readXML(p2, input2);

    EXPECT_EQ(648, simulation1->getContainer()->getParticleNumber());
    EXPECT_EQ(648, simulation2->getContainer()->getParticleNumber());

    simulation1->getContainer()->updateF(0);
    simulation2->getContainer()->updateF(1);

    for (int i = 0; i < 648; ++i) {
        EXPECT_NEAR(0, ArrayUtils::L2Norm(simulation1->getContainer()->getParticles()[i].getF() - simulation2->getContainer()->getParticles()[i].getF()), 1e-10);
    }

    simulation1->run();
    simulation2->run();

    for (int i = 0; i < 648; ++i) {
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ref[i].getX() - simulation1->getContainer()->getParticles()[i].getX()), 1e-10);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ref[i].getV() - simulation1->getContainer()->getParticles()[i].getV()), 1e-10);
        EXPECT_NEAR(ref[i].getM(), simulation1->getContainer()->getParticles()[i].getM(), 1e-10);
        EXPECT_NEAR(ref[i].getEpsilon(), simulation1->getContainer()->getParticles()[i].getEpsilon(), 1e-10);
        EXPECT_NEAR(ref[i].getSigma(), simulation1->getContainer()->getParticles()[i].getSigma(), 1e-10);
    }


    for (int i = 0; i < 648; ++i) {
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ref[i].getX() - simulation2->getContainer()->getParticles()[i].getX()), 1e-10);
        EXPECT_NEAR(0, ArrayUtils::L2Norm(ref[i].getV() - simulation2->getContainer()->getParticles()[i].getV()), 1e-10);
        EXPECT_NEAR(ref[i].getM(), simulation2->getContainer()->getParticles()[i].getM(), 1e-10);
        EXPECT_NEAR(ref[i].getEpsilon(), simulation2->getContainer()->getParticles()[i].getEpsilon(), 1e-10);
        EXPECT_NEAR(ref[i].getSigma(), simulation2->getContainer()->getParticles()[i].getSigma(), 1e-10);
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("Parallelization - Simulation test failed\n\n");
    } else {
        test_logger->info("Parallelization - Simulation test passed\n\n");
    }
}