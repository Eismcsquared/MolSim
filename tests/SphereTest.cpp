#include <spdlog/spdlog.h>
#include <gtest/gtest.h>
#include <vector>
#include "body/Sphere.h"
#include "Logger.h"
#include "utils/ArrayUtils.h"

class SphereTest: public ::testing::Test {
protected:
    Sphere sphere2D = Sphere({0, 0, 0}, {0, 0, 0}, 1000, 1, 1, 1, 2);
    Sphere sphere3D = Sphere({5, 5, 5}, {-1, 0, 1}, 100, 1, 3, 0.5, 3);
    const double pi = 3.141592653589793;
};

// Test whether a 2D sphere generates the correct particles in terms their position, mass and velocity including Brownian motion.
TEST_F(SphereTest, CreateParticles2D) {
    test_logger->info("Sphere - 2D Sphere test");
    std::vector<Particle> particles;
    sphere2D.createParticles(particles);
    double mean = 0;
    double meanSquire = 0;
    for(Particle &p: particles) {
        EXPECT_TRUE(ArrayUtils::L2Norm(p.getX()) <= 1000) << "Wrong particle created at position ("
        << p.getX()[0] << ", " << p.getX()[1] << ", " << p.getX()[2] << ")";
        EXPECT_NEAR(1, p.getM(), 1e-12) << "Particle with wrong mass encountered. Expected: 1, but got: " << p.getM();
        mean += ArrayUtils::L2Norm(p.getV());
        meanSquire += pow(ArrayUtils::L2Norm(p.getV()), 2);
    }
    mean /= particles.size();
    meanSquire /= particles.size();
    EXPECT_NEAR(pow(pi / 2, 0.5), mean, 2e-3) << "Wrong mean velocity of Brownian motion. Expected: "
    << pow(pi / 2, 0.5) << ", but got " << mean;
    EXPECT_NEAR(2, meanSquire, 5e-3) << "Wrong mean squared velocity of Brownian motion. Expected: "
    << 2 << ", but got " << meanSquire;
    if(::testing::Test::HasFailure()) {
        test_logger->info("Sphere - 2D Sphere test failed");
    } else {
        test_logger->info("Sphere - 2D Sphere test passed");
    }
}

// Test whether a 2D sphere generates the correct particles in terms their position, mass and velocity including Brownian motion.
TEST_F(SphereTest, CreateParticles3D) {
    test_logger->info("Sphere - 3D Sphere test");
    std::vector<Particle> particles;
    sphere3D.createParticles(particles);
    double mean = 0;
    double meanSquire = 0;
    for(Particle &p: particles) {
        EXPECT_TRUE(ArrayUtils::L2Norm(p.getX() - std::array<double, 3>{5, 5, 5}) <= 300) << "Wrong particle created at position ("
        << p.getX()[0] << ", " << p.getX()[1] << ", " << p.getX()[2] << ")";
        EXPECT_NEAR(1, p.getM(), 1e-12) << "Particle with wrong mass encountered. Expected: 1, but got: " << p.getM();
        mean += ArrayUtils::L2Norm(p.getV() - std::array<double, 3>{-1, 0, 1});
        meanSquire += pow(ArrayUtils::L2Norm(p.getV() - std::array<double, 3>{-1, 0, 1}), 2);
    }
    mean /= particles.size();
    meanSquire /= particles.size();
    EXPECT_NEAR(pow(2 / pi, 0.5), mean, 1e-3) << "Wrong mean velocity of Brownian motion. Expected: "
    << pow(2 / pi, 0.5) << ", but got " << mean;
    EXPECT_NEAR(0.75, meanSquire, 2.5e-3) << "Wrong mean squared velocity of Brownian motion. Expected: "
    << 0.75 << ", but got " << meanSquire;
    if(::testing::Test::HasFailure()) {
        test_logger->info("Sphere - 3D Sphere test failed");
    } else {
        test_logger->info("Sphere - 3D Sphere test passed");
    }
}