#include "Container.h"

Container::Container(std::vector<Particle> &particles, std::unique_ptr<Force>& f_ptr): particles(particles), f(*f_ptr){}

unsigned long Container::getParticleNumber() const {
    return particles.size();
}

std::vector<Particle> &Container::getParticles() const {
    return particles;
}

void Container::addParticle(const Particle &particle) {
    this->particles.push_back(particle);
    this->Container::updateF();
}

void Container::addCuboid(const Cuboid &cuboid) {
    cuboid.createParticles(particles);
    this->Container::updateF();
}

void Container::updateF() {
    updateF(true);
}

void Container::simulate(double end_time, double delta_t, const std::string &out_name, const std::string &output_format,
                         unsigned int output_frequency, bool save_output, bool newton3) {
    int max_iteration = (int) (end_time / delta_t);
    for (int iteration = 0; iteration < max_iteration; iteration++) {
        if (iteration % output_frequency == 0 && save_output) {
            plotParticles(iteration, out_name, output_format);
        }
        // Calculate the position
        updateX(delta_t);
        // Calculate the force
        updateF(newton3);
        // Calculate the velocity
        //updateV(delta_t);
        // Plot every 10th iteration
        spdlog::trace("Iteration {} finished.", iteration + 1);
    }
    spdlog::trace("output written. Terminating...");
}

void Container::plotParticles(int iteration, const std::string &out_name, const std::string &output_format) {
    spdlog::trace("Plotting Particles...");

    if(output_format == "vtu") {
        outputWriter::VTKWriter writer;
        writer.writeFile(out_name, iteration,particles);
    }
    else if(output_format == "xyz") {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particles, out_name, iteration);
    }
}


std::string Container::toString() {
    std::stringstream buf;
    buf << "Number of particles: " << particles.size() << std::endl;
    for (auto &p : particles) {
        buf << p.toString() << std::endl;
    }
    return buf.str();
}

bool Container::operator==(const Container &other) const {
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

std::unique_ptr<Iterator> Container::iterator() const {
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


