#include <gtest/gtest.h>
#include <vector>
#include <spdlog/spdlog.h>
#include "container/LinkedCellContainer.h"
#include "force/LennardJonesForce.h"

class ThermostatTest: public ::testing::Test {

protected:
    std::vector<Particle> particles1;
    std::vector<Particle> particles2;
    std::shared_ptr<spdlog::logger> test_logger = spdlog::get("test_logger");

    void SetUp() override {

        double sqrt3over2 = sqrt(3) / 2;

        // A small 2D example with few particles. Initial temperature is 21.
        particles1.emplace_back(std::array<double, 3>{0, 0, 0}, std::array<double, 3>{10, 10, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{10, 0, 0}, std::array<double, 3>{17, 10, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{5, 10 * sqrt3over2, 0}, std::array<double, 3>{13.5, 10 + sqrt3over2 * 7, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{-5, 10 * sqrt3over2, 0}, std::array<double, 3>{6.5, 10 + sqrt3over2 * 7, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{-10, 0, 0}, std::array<double, 3>{3, 10, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{-5, -10 * sqrt3over2, 0}, std::array<double, 3>{6.5, 10 - sqrt3over2 * 7, 0}, 1);
        particles1.emplace_back(std::array<double, 3>{5, -10 * sqrt3over2, 0}, std::array<double, 3>{13.5, 10 - sqrt3over2 * 7, 0}, 1);

        // A large 3D example with 200000 particles. Initial temperature is 25.
        Cuboid cuboid1({10, 10, 10}, {0, 0, 0}, {100, 100, 10}, 1, pow(2, 1.0 / 6), 5, 3);
        Cuboid cuboid2({10, 150, 10}, {0, 0, 0}, {100, 100, 10}, 2, pow(2, 1.0 / 6), 5 / sqrt(2), 3);
        cuboid1.createParticles(particles2);
        cuboid2.createParticles(particles2);
    }
};

// Test whether the helper function temperature computes the temperature correctly, which will be used in the tests.
TEST_F(ThermostatTest, Temperature) {
    test_logger->info("Thermostat - Temperature computation test");

    Thermostat thermostat2D(10, 1, 10, 2);
    Thermostat thermostat3D(10, 1, 10, 3);

    EXPECT_FLOAT_EQ(21, Thermostat::temperature(particles1, 2));
    EXPECT_NEAR(25, Thermostat::temperature(particles2, 3), .25);

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Temperature computation test failed\n\n");
    } else {
        test_logger->info("Thermostat - Temperature computation test passed\n\n");
    }
}

// Test the heating process by the thermostats, both cases maxDelta > target_T - initial_T and maxDelta < target_T - initial_T are considered.
TEST_F(ThermostatTest, Heating) {

    test_logger->info("Thermostat - Heating test");

    Thermostat heater2D(84, 1, std::numeric_limits<double>::infinity(), 2);
    Thermostat gradualHeater2D(30, 1, 7, 2);

    gradualHeater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{10 + 14 / sqrt(3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{10 + 7 / sqrt(3), 17, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{10 - 7 / sqrt(3), 17, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{10 - 14 / sqrt(3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{10 - 7 / sqrt(3), 3, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{10 + 7 / sqrt(3), 3, 0}), 1e-12);

    heater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{24, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{17, 10 + 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{3, 10 + 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{-4, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{3, 10 - 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{17, 10 - 7 * sqrt(3), 0}), 1e-12);

    Thermostat heater3D(30, 1, 10, 3);
    Thermostat gradualHeater3D(50, 1, 10, 3);

    heater3D.apply(particles2);

    EXPECT_FLOAT_EQ(30, Thermostat::temperature(particles2, 3));

    gradualHeater3D.apply(particles2);

    EXPECT_FLOAT_EQ(40, Thermostat::temperature(particles2, 3));

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Heating test failed\n\n");
    } else {
        test_logger->info("Thermostat - Heating test passed\n\n");
    }
}

// Test the cooling process by the thermostats, both cases maxDelta > initial_T - target_T and maxDelta < initial_T - target_T are considered.
TEST_F(ThermostatTest, Cooling) {
    test_logger->info("Thermostat - Cooling test");

    Thermostat heater2D(7, 1, std::numeric_limits<double>::infinity(), 2);
    Thermostat gradualHeater2D(10, 1, 7, 2);

    gradualHeater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{10 + 7 * sqrt(2.0 / 3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{10 + 3.5 * sqrt(2.0 / 3), 10 + 3.5 * sqrt(2), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{10 - 3.5 * sqrt(2.0 / 3), 10 + 3.5 * sqrt(2), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{10 - 7 * sqrt(2.0 / 3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{10 - 3.5 * sqrt(2.0 / 3), 10 - 3.5 * sqrt(2), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{10 + 3.5 * sqrt(2.0 / 3), 10 - 3.5 * sqrt(2), 0}), 1e-12);

    heater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{10 + 7 / sqrt(3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{10 + 3.5 / sqrt(3), 13.5, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{10 - 3.5 / sqrt(3), 13.5, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{10 - 7 / sqrt(3), 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{10 - 3.5 / sqrt(3), 6.5, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{10 + 3.5 / sqrt(3), 6.5, 0}), 1e-12);

    Thermostat heater3D(20, 1, 10, 3);
    Thermostat gradualHeater3D(10, 1, 5, 3);

    heater3D.apply(particles2);

    EXPECT_FLOAT_EQ(20, Thermostat::temperature(particles2, 3));

    gradualHeater3D.apply(particles2);

    EXPECT_FLOAT_EQ(15, Thermostat::temperature(particles2, 3));

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Cooling test failed\n\n");
    } else {
        test_logger->info("Thermostat - Cooling test passed\n\n");
    }
}

// Test the thermostat holding temperature constant.
TEST_F(ThermostatTest, HoldingConstant) {
    test_logger->info("Thermostat - Holding temperature constant test");

    Thermostat heater2D(21, 1, std::numeric_limits<double>::infinity(), 2);

    heater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{17, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{13.5, 10 + 3.5 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{6.5, 10 + 3.5 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{3, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{6.5, 10 - 3.5 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{13.5, 10 - 3.5 * sqrt(3), 0}), 1e-12);

    Thermostat heater3D(25, 1, 1, 3);

    heater3D.apply(particles2);

    EXPECT_FLOAT_EQ(25, Thermostat::temperature(particles2, 3));

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Holding temperature constant test failed\n\n");
    } else {
        test_logger->info("Thermostat - Holding temperature constant test passed\n\n");
    }
}

// Add some stationary particles and test whether the effect of the thermostat stays the same.
TEST_F(ThermostatTest, Stationary) {
    test_logger->info("Thermostat - Including stationary particles test");

    double sqrt3over2 = sqrt(3) / 2;
    particles1.emplace_back(std::array<double, 3>{5, 0, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
    particles1.emplace_back(std::array<double, 3>{2.5, 5 * sqrt3over2, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
    particles1.emplace_back(std::array<double, 3>{-2.5, 5 * sqrt3over2, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
    particles1.emplace_back(std::array<double, 3>{-5, 0, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
    particles1.emplace_back(std::array<double, 3>{-2.5, -5 * sqrt3over2, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);
    particles1.emplace_back(std::array<double, 3>{2.5, -5 * sqrt3over2, 0}, std::array<double, 3>{0, 0, 0}, 1, 0, 5, 1, true);

    Thermostat heater2D(84, 1, std::numeric_limits<double>::infinity(), 2);

    EXPECT_FLOAT_EQ(21, Thermostat::temperature(particles1, 2));

    heater2D.apply(particles1);

    EXPECT_LE(ArrayUtils::L2Norm(particles1[0].getV() - std::array<double, 3>{10, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[1].getV() - std::array<double, 3>{24, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[2].getV() - std::array<double, 3>{17, 10 + 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[3].getV() - std::array<double, 3>{3, 10 + 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[4].getV() - std::array<double, 3>{-4, 10, 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[5].getV() - std::array<double, 3>{3, 10 - 7 * sqrt(3), 0}), 1e-12);
    EXPECT_LE(ArrayUtils::L2Norm(particles1[6].getV() - std::array<double, 3>{17, 10 - 7 * sqrt(3), 0}), 1e-12);

    for (int i = 7; i < 13; ++i) {
        EXPECT_FLOAT_EQ(0, ArrayUtils::L2Norm(particles1[i].getV()));
    }

    Cuboid cuboid3({10, 290, 10}, {0, 0, 0}, {100, 100, 10}, 1, pow(2, 1.0 / 6), 0, 3, 0, 5, 1, true);
    cuboid3.createParticles(particles2);

    EXPECT_NEAR(25, Thermostat::temperature(particles2, 3), .25);

    Thermostat heater3D(30, 1, 10, 3);

    heater3D.apply(particles2);

    EXPECT_FLOAT_EQ(30, Thermostat::temperature(particles2, 3));

    for (int i = 200000; i < 300000; ++i) {
        EXPECT_FLOAT_EQ(0, ArrayUtils::L2Norm(particles2[i].getV()));
    }

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Including stationary particles test failed\n\n");
    } else {
        test_logger->info("Thermostat - Including stationary particles test passed\n\n");
    }
}

// Test the integration of thermostats into particle containers.
TEST_F(ThermostatTest, IntegrationContainer) {
    test_logger->info("Thermostat - Holding temperature constant test");

    std::vector<Particle> particles;
    Cuboid cuboid1({10, 10, 10}, {0, 0, 0}, {10, 10, 10}, 1, pow(2, 1.0 / 6), 5, 3);
    Cuboid cuboid2({10, 30, 10}, {0, 0, 0}, {10, 10, 10}, 2, pow(2, 1.0 / 6), 5 / sqrt(2), 3);
    cuboid1.createParticles(particles);
    cuboid2.createParticles(particles);
    std::unique_ptr<Force> force = std::make_unique<LennardJonesForce>();
    LinkedCellContainer container(particles, force, {30, 60, 30}, 3, {REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING, REFLECTING});

    std::unique_ptr<Thermostat> holdingConstant = std::make_unique<Thermostat>(25, 5, 10, 3);
    std::unique_ptr<Thermostat> heating = std::make_unique<Thermostat>(50, 5, 3, 3);
    std::unique_ptr<Thermostat> cooling = std::make_unique<Thermostat>(10, 5, 7, 3);

    container.setThermostat(holdingConstant);
    container.simulate(5e-5, 1e-5, "", "vtu", 10, false);
    EXPECT_FLOAT_EQ(25, Thermostat::temperature(particles, 3));

    // gradual heating
    container.setThermostat(heating);
    container.simulate(2.5e-4, 1e-5, "", "vtu", 10, false);
    EXPECT_NEAR(37, Thermostat::temperature(particles, 3), 1e-2);
    container.simulate(6.5e-4, 1e-5, "", "vtu", 10, false);
    EXPECT_FLOAT_EQ(50, Thermostat::temperature(particles, 3));

    // gradual cooling
    container.setThermostat(cooling);
    container.simulate(8.5e-4, 1e-5, "", "vtu", 10, false);
    EXPECT_NEAR(22, Thermostat::temperature(particles, 3), 1e-2);
    container.simulate(1.25e-3, 1e-5, "", "vtu", 10, false);
    EXPECT_FLOAT_EQ(10, Thermostat::temperature(particles, 3));

    if (::testing::Test::HasFailure()) {
        test_logger->info("Thermostat - Holding temperature constant test failed\n\n");
    } else {
        test_logger->info("Thermostat - Holding temperature constant test passed\n\n");
    }

}
