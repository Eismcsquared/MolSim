#include <utility>
#include "Simulation.h"
#include "outputWriter/StateWriter.h"



Simulation::Simulation(std::unique_ptr<ParticleContainer> &container, double endTime, double deltaT,
                       std::string outputFile, std::string outputFormat, unsigned int outputFrequency)
        : container(std::move(container)), endTime(endTime), deltaT(deltaT), outputFormat(std::move(outputFormat)), outputFile(std::move(outputFile)),
          outputFrequency(outputFrequency), saveOutput(true), statistics(nullptr) {}

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

void Simulation::setCheckpointingFile(const std::string &checkpointingFile) {
    Simulation::checkpointingFile = checkpointingFile;
}

void Simulation::setStatistics(std::shared_ptr<Statistics> statistics) {
    Simulation::statistics = statistics;
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

std::string Simulation::getCheckpointingFile() const {
    return checkpointingFile;
}

Statistics &Simulation::getStatistics() const {
    return *statistics;
}

void Simulation::run() {
    container->simulate(endTime, deltaT, outputFile, outputFormat, outputFrequency, saveOutput, statistics);

    if (!checkpointingFile.empty()) {
        StateWriter::saveState(container->getParticles(), checkpointingFile);
    }
}

