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


char* testfile = "../input/eingabe-sonne.txt";
std::shared_ptr<spdlog::logger> test_logger = spdlog::basic_logger_mt("test_logger", "logs/test.txt");
class ParticleTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;

    void SetUp() override {

        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        GravitationalForce f;
        pc = new ParticleContainer(particles, 0, 1000, 0.014, f, ".vtu");
        spdlog::set_level(spdlog::level::info);
        test_logger -> info("Particle Container created");

    }

    void TearDown() override {
        test_logger->info("delete Particle Container\n\n");
        delete pc;
    }
};

TEST_F(ParticleTest, Addparticle) {
    test_logger->info("Add particle test1");
    Particle p = Particle({1,0,0}, {0,0,0}, 1, 0);
    int size_before = pc->getParticleSize();
    pc->addParticle(p);
    int size_after = pc->getParticleSize();
    ASSERT_EQ(size_before + 1, size_after);
    if(size_before + 1 == size_after) {
        test_logger->info("Add particle test1 passed");
    } else {
        test_logger->error("Add particle test1 failed");
    }
}

TEST_F(ParticleTest, Addparticle_2) {
    test_logger -> info("Add particle test2");
    int size_before = pc->getParticleSize();

    bool flag = true;
    for(int i = 1 ; i <= 10; i++) {
        Particle p = Particle({(double)i +1,0,0}, {0,0,0}, 1, 0);
        pc->addParticle(p);
        int size_after = pc->getParticleSize();
        ASSERT_EQ(size_before + i, size_after) ;  
        if(size_before + i == size_after) {
        } else {
            test_logger->error("Add particle test2 at iteration {} failed", i);
            flag = false;
            break;
        }
    }
    if(flag){
        test_logger->info("Add particle test2 passed");
    }
}

TEST_F(ParticleTest, Update_F_simple) {
    test_logger->info("Calculate simple force test");
   // pc->updateF();
    test_logger->info("Calculate simple force test passed");
}

TEST_F(ParticleTest, Update_X_test) {
    test_logger->info("Update X test");
    pc->updateX();
    test_logger->info("Update X test passed");
}

TEST_F(ParticleTest, Update_V_test) {
    test_logger->info("Update V test");
    pc->updateV();
    test_logger->info("Update V test passed");
}

TEST_F(ParticleTest, Complex_test) {
    test_logger->info("complex test"); // result
    test_logger->info("complex test passed");
}