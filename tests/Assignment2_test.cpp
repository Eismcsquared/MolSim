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

TEST(Test3, Test3) {
    std::cout << "Test1" << std::endl;
    ASSERT_TRUE(true);
}

TEST(Test4, Test4) {
    std::cout << "Test2" << std::endl;
    EXPECT_FALSE(false);
}