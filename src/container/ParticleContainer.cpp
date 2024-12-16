#include "ParticleContainer.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/StateWriter.h"

ParticleContainer::ParticleContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr) :
        ParticleContainer(particles, f_ptr, 0) {}

ParticleContainer::ParticleContainer(std::vector<Particle> &particles,
                                             std::unique_ptr<Force> &f_ptr, double g) :
        particles(particles), force(std::move(f_ptr)), g(g), particleNumber(particles.size()) {

    for (Particle &p: particles) {
        if (!p.isInDomain()) {
            particleNumber--;
        }
    }
}


unsigned long ParticleContainer::getParticleNumber() const {
    return particleNumber;
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
    if (particle.isInDomain()) {
        this->particles.push_back(particle);
        particleNumber++;
    }
}

void ParticleContainer::addCluster(const Cluster &cluster) {
    unsigned long sizeOld = particles.size();
    cluster.createParticles(particles);
    particleNumber += particles.size() - sizeOld;
}

void ParticleContainer::simulate(double start_time, double end_time, double delta_t, const std::string &out_name, const std::string &output_format,
                                 unsigned int output_frequency, bool save_output) {
    int start_iteration = static_cast<int>(std::round(start_time / delta_t));
    int end_iteration = static_cast<int>(std::round(end_time / delta_t));

    // compute initial forces.
    updateF();
    // save the initial state also.
    if (save_output && start_iteration % output_frequency == 0) {
        plotParticles(start_iteration, out_name, output_format);
    }

    long double duration = 0;
    unsigned long numberUpdates = 0;

    for (int iteration = start_iteration + 1; iteration <= end_iteration; iteration++) {

        // Start the timer.
        auto start = std::chrono::high_resolution_clock::now();

        // Calculate the position
        updateX(delta_t);

        // Calculate the force
        updateF();

        // Calculate the velocity
        updateV(delta_t);

        if (thermostat && iteration % thermostat->getPeriode() == 0) {
            thermostat->apply(particles);
        }

        // Stop the timer
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        duration += (end - start).count();
        numberUpdates += particleNumber;

        if (iteration % output_frequency == 0 && save_output) {
            plotParticles(iteration, out_name, output_format);
        }

        spdlog::trace("Iteration {} finished.", iteration);
    }

    spdlog::trace("output written. Terminating...");

    // print measurement results.
    spdlog::set_level(spdlog::level::info);
    spdlog::info("Duration: {}s", duration);
    spdlog::info("Updates per second: {}", numberUpdates / duration);
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


