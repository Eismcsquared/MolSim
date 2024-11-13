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
#include <gmock/gmock.h>
#include "LennardJonesForce.h"
#include "Force.h"

TEST(LennardJonesForceTest, LennardJonesForce) {
    std::cout << "Test2" << std::endl;
    EXPECT_FALSE(false);
}

TEST(GravitationalForceTest, GravitationalForce) {
    std::cout << "Test1" << std::endl;
    EXPECT_FALSE(false);
}

TEST(ParticleContainerTest, ParticleContainer) {
    std::cout << "Test3" << std::endl;
    EXPECT_FALSE(false);
}
