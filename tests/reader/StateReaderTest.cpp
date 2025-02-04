#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <filesystem>
#include "outputWriter/StateWriter.h"
#include "inputReader/StateReader.h"
#include "body/Cuboid.h"

class StateReaderTest: public ::testing::Test {
protected:
    std::string fileName = "../tests/test_cases/phase_space.txt";
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");
    std::vector<Particle> particles;

    void SetUp() override {
        Cuboid c1(std::array<double, 3>{0, 0, 0}, std::array<double, 3>{1, 0, 5}, std::array<unsigned int, 3>{10, 10, 10}, 1, 1.5, 10, 3, 0, 5, 1);
        Cuboid c2(std::array<double, 3>{15, 20, -15}, std::array<double, 3>{6, -3, 2}, std::array<unsigned int, 3>{10, 10, 10}, 2, 1.2, 10, 3, 1, 3, 1.2);
        c1.createParticles(particles);
        c2.createParticles(particles);
    }

};


// Test whether the state reader loads the state of particles correctly.
TEST_F(StateReaderTest, ReloadState) {
    test_logger->info("StateReaderTest - Reload states test");

    std::filesystem::remove(fileName);
    StateWriter::saveState(particles, fileName);
    std::vector<Particle> loadedParticles;
    StateReader::loadState(loadedParticles, fileName);

    EXPECT_EQ(2000, loadedParticles.size());
    for (int i = 0; i < 2000; ++i) {
        EXPECT_EQ(particles[i], loadedParticles[i]);
    }
    std::filesystem::remove(fileName);
    if (::testing::Test::HasFailure()) {
        test_logger->info("StateReaderTest - Reload states test failed\n\n");
    } else {
        test_logger->info("StateReaderTest - Reload states test passed\n\n");
    }
}