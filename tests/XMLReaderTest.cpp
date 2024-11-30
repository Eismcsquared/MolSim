#include <gtest/gtest.h>
#include <string>
#include "inputReader/XMLReader.h"
#include "Logger.h"

class XMLReaderTest: public ::testing::Test {
protected:
    std::string inputFile;
    std::unique_ptr<Simulation> simulation;
};

TEST_F(XMLReaderTest, Assigment1Input) {
    inputFile = "../tests/test_cases/assignment1.xml";
    simulation = XMLReader::readXML(inputFile);
    std::cout << simulation->getContainer()->toString() << std::endl;
    ASSERT_EQ(4, simulation->getContainer()->getParticleNumber());
}
