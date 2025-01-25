#include <vector>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include "body/Membrane.h"

class MembraneTest: public ::testing::Test {
protected:
    Membrane membrane1 = Membrane(std::array<double, 3>{1, 3, 0}, std::array<double, 3>{1, 2, 1}, std::array<unsigned int, 2>{10, 20},
                                 1, 1, 0, 2, 0, 1, 1, 100, 1);
    Membrane membrane2 = Membrane(std::array<double, 3>{7, 7, 5}, std::array<double, 3>{0, 0, 0}, std::array<unsigned int, 2>{20, 10},
                                  1, 1, 0, 2, 0, 1, 1, 100, 1);
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");
};

// Test whether the membrane generates particles with the given parameters and the correct neighbourhood relationship.
TEST_F(MembraneTest, CreateParticle) {
    test_logger->info("Membrane - Create particles test");
    std::vector<Particle> particles;
    membrane1.createParticles(particles);
    EXPECT_EQ(200, particles.size());
    for (Particle &particle1 : particles) {
        EXPECT_FLOAT_EQ(0, ArrayUtils::L2NormSquare(particle1.getV() - std::array<double, 3>{1, 2, 1}));
        EXPECT_FLOAT_EQ(1, particle1.getM());
        EXPECT_EQ(0, particle1.getType());
        EXPECT_FLOAT_EQ(1, particle1.getEpsilon());
        EXPECT_FLOAT_EQ(1, particle1.getSigma());
        EXPECT_FLOAT_EQ(100, particle1.getK());
        EXPECT_FLOAT_EQ(1, particle1.getR0());

        // particles are neighbours if their distance is 1 and diagonal neighbours if their distance is sqrt(2).
        for (Particle &particle2 : particles) {
            if (std::abs(ArrayUtils::L2NormSquare(particle1.getX() - particle2.getX()) - 1) < 1e-12) {
                EXPECT_TRUE(particle1.isNeighbour(particle2));
            } else {
                EXPECT_FALSE(particle1.isNeighbour(particle2));
            }
            if (std::abs(ArrayUtils::L2NormSquare(particle1.getX() - particle2.getX()) - 2) < 1e-12) {
                EXPECT_TRUE(particle1.isDiagonalNeighbour(particle2));
            } else {
                EXPECT_FALSE(particle1.isDiagonalNeighbour(particle2));
            }
        }
    }
    membrane2.createParticles(particles);
    EXPECT_EQ(400, particles.size());

    for (int i = 200; i < 400; ++i) {
        // particles from different membranes are never neighbours.
        for (int j = 0; j < 200; ++j) {
            EXPECT_FALSE(particles[i].isNeighbour(particles[j]));
            EXPECT_FALSE(particles[j].isNeighbour(particles[i]));
            EXPECT_FALSE(particles[i].isDiagonalNeighbour(particles[j]));
            EXPECT_FALSE(particles[j].isDiagonalNeighbour(particles[i]));
        }

        for (int j = 200; j < 400; ++j) {
            if (std::abs(ArrayUtils::L2NormSquare(particles[i].getX() - particles[j].getX()) - 1) < 1e-12) {
                EXPECT_TRUE(particles[i].isNeighbour(particles[j]));
            } else {
                EXPECT_FALSE(particles[i].isNeighbour(particles[j]));
            }
            if (std::abs(ArrayUtils::L2NormSquare(particles[i].getX() - particles[j].getX()) - 2) < 1e-12) {
                EXPECT_TRUE(particles[i].isDiagonalNeighbour(particles[j]));
            } else {
                EXPECT_FALSE(particles[i].isDiagonalNeighbour(particles[j]));
            }
        }
    }
    if (::testing::Test::HasFailure()) {
        test_logger->info("Membrane - Create particles test failed\n\n");
    } else {
        test_logger->info("Membrane - Create particles test passed\n\n");
    }
}
