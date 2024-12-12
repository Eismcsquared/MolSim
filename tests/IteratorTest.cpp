#include <gtest/gtest.h>
#include <memory>
#include "Logger.h"
#include "inputReader/FileReader.h"
#include "body/Particle.h"
#include "container/DirectSumContainer.h"
#include "force/LennardJonesForce.h"

class IteratorTest: public ::testing::Test {
protected:
    std::unique_ptr<std::vector<Particle>> particles;
    std::unique_ptr<DirectSumContainer> pc;
    std::unique_ptr<Force> f;
    char* testfile = const_cast<char*>("../tests/test_cases/assignment1.txt");


    void SetUp() override {
        FileReader fileReader;
        particles = std::make_unique<std::vector<Particle>>();
        fileReader.readFile(*particles, testfile);

        f = std::make_unique<LennardJonesForce>();
        pc = std::make_unique<DirectSumContainer>(particles, f);
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("ParticleContainer created");
    }

    void TearDown() override {
        test_logger->info("ParticleContainer deleted\n\n");
    }
};

//Test whether the iterator of the particle container works correctly
TEST_F(IteratorTest, ParticleContainerIter) {
    test_logger->info("Particle container iterator test");
    std::unique_ptr<Iterator> cur = pc->iterator();
    for (int i = 0; i < pc->getParticleNumber(); ++i) {
        Particle current = cur->next();
        if (!(pc->getParticles()[i] == current)) {
            test_logger->error("Particle container iterator test failed");
        }
        EXPECT_EQ(pc->getParticles()[i], current);
    }
    if (cur->hasNext()) {
        test_logger->error("Particle container iterator test failed");
    }
    ASSERT_FALSE(cur->hasNext());
    test_logger->info("Particle container iterator test passed");
}
