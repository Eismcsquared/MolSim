#include <utility>
#include "Simulation.h"


Simulation::Simulation(std::unique_ptr<ParticleContainer> &container, double endTime, double deltaT,
                       std::string outputFile, std::string outputFormat, unsigned int outputFrequency)
        : container(std::move(container)), endTime(endTime), deltaT(deltaT), outputFormat(std::move(outputFormat)), outputFile(std::move(outputFile)),
          outputFrequency(outputFrequency), saveOutput(true), newton3(true) {}


void Simulation::setEndTime(double endTime) {
    Simulation::endTime = endTime;
}

void Simulation::setDeltaT(double deltaT) {
    Simulation::deltaT = deltaT;
}

void Simulation::setOutputFormat(const std::string &outputFormat) {
    Simulation::outputFormat = outputFormat;
}

void Simulation::setOutputFile(const std::string &outputFile) {
    Simulation::outputFile = outputFile;
}

void Simulation::setOutputFrequency(unsigned int outputFrequency) {
    Simulation::outputFrequency = outputFrequency;
}

void Simulation::setSaveOutput(bool saveOutput) {
    Simulation::saveOutput = saveOutput;
}

void Simulation::setNewton3(bool newton3) {
    Simulation::newton3 = newton3;
}


const std::unique_ptr<ParticleContainer> &Simulation::getContainer() const {
    return container;
}

double Simulation::getEndTime() const {
    return endTime;
}

double Simulation::getDeltaT() const {
    return deltaT;
}

const std::string &Simulation::getOutputFile() const {
    return outputFile;
}

const std::string &Simulation::getOutputFormat() const {
    return outputFormat;
}

unsigned int Simulation::getOutputFrequency() const {
    return outputFrequency;
}

bool Simulation::isSaveOutput() const {
    return saveOutput;
}

bool Simulation::isNewton3() const {
    return newton3;
}


void Simulation::run() {
    container->simulate(endTime, deltaT, outputFile, outputFormat, outputFrequency, saveOutput, newton3);
}
