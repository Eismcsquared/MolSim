#include "ParticleContainer.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/StateWriter.h"

ParticleContainer::ParticleContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr) :
        particles(particles), force(std::move(f_ptr)), g(0) {}

ParticleContainer::ParticleContainer(std::vector<Particle> &particles,
                                             std::unique_ptr<Force> &f_ptr, double g) :
        particles(particles), force(std::move(f_ptr)), g(g) {}


unsigned long ParticleContainer::getParticleNumber() const {
    int num = 0;
    for (const auto& p: getParticles()) {
        if (p.isInDomain()) {
            num++;
        }
    }
    return num;
}

std::vector<Particle> &ParticleContainer::getParticles() const {
    return particles;
}

double ParticleContainer::getG() const {
    return g;
}

std::unique_ptr<Thermostat> &ParticleContainer::getThermostat(){
    return thermostat;
}

void ParticleContainer::setF(std::unique_ptr<Force> &f) {
    this->force = std::move(f);
}


void ParticleContainer::setG(double g) {
    ParticleContainer::g = g;
}

void ParticleContainer::setThermostat(std::unique_ptr<Thermostat> &thermostat) {
    ParticleContainer::thermostat = std::move(thermostat);
}

void ParticleContainer::addParticle(const Particle &particle) {
    this->particles.push_back(particle);
    this->ParticleContainer::updateF();
}

void ParticleContainer::addCluster(const Cluster &cluster) {
    cluster.createParticles(particles);
    this->ParticleContainer::updateF();
}

void ParticleContainer::updateF() {
    updateF(true);
}

void ParticleContainer::simulate(double end_time, double delta_t, const std::string &out_name, const std::string &output_format,
                                 unsigned int output_frequency, bool save_output, const std::string& checkpointing, bool newton3) {
    int max_iteration = static_cast<int>(std::round(end_time / delta_t));

    // save the initial state also.
    if (save_output) {
        plotParticles(0, out_name, output_format);
    }

    for (int iteration = 1; iteration <= max_iteration; iteration++) {

        // Calculate the position
        updateX(delta_t);

        // Calculate the force
        updateF(newton3);

        // Calculate the velocity
        updateV(delta_t);

        if (iteration % output_frequency == 0 && save_output) {
            plotParticles(iteration, out_name, output_format);
        }

        if (thermostat && iteration % thermostat->getPeriode() == 0) {
            thermostat->apply(particles);
        }

        spdlog::trace("Iteration {} finished.", iteration);
    }

    if (!checkpointing.empty()) {
        StateWriter::saveState(particles, checkpointing);
    }

    spdlog::trace("output written. Terminating...");
}

void ParticleContainer::plotParticles(int iteration, const std::string &out_name, const std::string &output_format) {
    spdlog::trace("Plotting Particles...");

    if(output_format == "vtu") {
        outputWriter::VTKWriter writer;
        writer.writeFile(out_name, iteration, particles);
    }
    else if(output_format == "xyz") {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particles, out_name, iteration);
    }
}


std::string ParticleContainer::toString() {
    std::stringstream buf;
    buf << "Number of particles: " << particles.size() << std::endl;
    for (auto &p : particles) {
        buf << p.toString() << std::endl;
    }
    return buf.str();
}

bool ParticleContainer::operator==(const ParticleContainer &other) const {
    if (getParticleNumber() != other.getParticleNumber()) {
        return false;
    }
    for (int i = 0; i < getParticleNumber(); ++i) {
        if (!(getParticles()[i] == other.getParticles()[i])) {
            return false;
        }
    }
    return true;
}

std::unique_ptr<Iterator> ParticleContainer::iterator() const {
    return std::make_unique<Iterator>(particles.begin(), particles.end());
}


Iterator::Iterator(std::vector<Particle>::iterator begin, std::vector<Particle>::iterator end): current(begin), end(end){}

Particle &Iterator::next() {
    Particle& p = *current;
    ++current;
    return p;
}

bool Iterator::hasNext() {
    return current != end;
}


