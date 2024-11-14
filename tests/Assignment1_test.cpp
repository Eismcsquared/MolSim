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
#include <fstream>
#include <sstream>




std::string TxttoString(const std::string& filePath) {
    std::ifstream file(filePath);
    std::stringstream buf;
    buf << file.rdbuf();
    return buf.str();
}

class ParticleTest : public ::testing::Test {
protected:
    std::vector<Particle> particles;
    ParticleContainer *pc;
    char* testfile = const_cast<char*>("../input/eingabe-sonne.txt");
    std::shared_ptr<spdlog::logger> test_logger = spdlog::basic_logger_mt("test_logger", "logs/test.txt");

    void SetUp() override {
        FileReader fileReader;
        fileReader.readFile(particles, testfile);
        pc = new ParticleContainer(particles, new GravitationalForce());
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
    bool flag = true;
    pc->updateF();
    pc->writeoutput("updateF_simple.txt");

    std::string actual = TxttoString("updateF_simple.txt");
    std::string expected = TxttoString("../tests/Answer_Ref/updateF_simple.txt");

    if(expected != actual) {
        test_logger->error("Calculate simple force test failed");
    }
    ASSERT_EQ(expected, actual) << "Calculate simple force test failed";
    if(flag){
        test_logger->info("Calculate simple force test passed");
    }
}

TEST_F(ParticleTest, Simulation_simple) {
    test_logger->info("Simulation simple test");
    bool flag = true;
    for(int iteration = 0 ; iteration <= 3; ++iteration) {
        pc->updateX(0.014);
        pc->updateF();
        pc->updateV(0.014);
    }
    std :: string actual = pc->writeoutput("simulation_simple.txt");
    std::string expected = TxttoString("../tests/Answer_Ref/Ans_simulation_simple.txt");

    if(expected != actual) {
        test_logger->error("Simulation simple test failed");
    }
    ASSERT_EQ(expected, actual) << "Simulation simple test failed";
    if(flag){
        test_logger->info("Simulation simple test passed");
    }
}


TEST_F(ParticleTest, Complex_simulation) {
    test_logger->info("Complex simulation test");
    bool flag = true;
    int max_iteration = (int) ((1000 - 0) / 0.014);
    for(int iteration = 0 ; iteration <= max_iteration; ++iteration) {
        pc->updateX(0.014);
        pc->updateF();
        pc->updateV(0.014);
    }
    std::string actual = pc->writeoutput("simulation_complex.txt");
    std::string expected = TxttoString("../tests/Answer_Ref/Ans_simulation_complex.txt");

    if(expected != actual) {
        test_logger->error("Complex simulation test failed");
    }
    ASSERT_EQ(expected, actual) << "Complex simulation test failed";
    if(flag){
        test_logger->info("Complex simulation test passed");
    }
}