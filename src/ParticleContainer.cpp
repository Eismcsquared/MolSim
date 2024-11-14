#include <iostream>
#include <vector>
#include <array>
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include <filesystem>


ParticleContainer::ParticleContainer(std::vector<Particle>& particles, Force* f)
    : particles(particles),
      f(*f) { // std::move if outputFormat is temporary
    // compute initial forces
    updateF();
    std::cout << "ParticleContainer generated!\n";
}



ParticleContainer::~ParticleContainer(){
    std::cout << "ParticleContainer destructed!\n";
    delete &f;
}


void ParticleContainer::addParticle(const Particle& particle){
    this->particles.push_back(particle);
    this->updateF();
}

void ParticleContainer::addCuboid(const Cuboid& cuboid) {
    cuboid.createParticles(particles);
    this->updateF();
}


void ParticleContainer::updateF(bool newton3) {
    if (newton3) {
        std::vector<std::array<double, 3>> newForces(particles.size(), {0, 0, 0});
        for (unsigned long i = 0; i < particles.size(); ++i) {
            for (unsigned long j = i + 1; j < particles.size(); ++j) {
                std::array<double, 3> forceIJ = f.force(particles[i], particles[j]);
                // update forces using the third Newton axiom
                newForces[i] = newForces[i] - forceIJ;
                newForces[j] = newForces[j] + forceIJ;
            }
        }
        for (unsigned long i = 0; i < particles.size(); ++i) {
            particles[i].setOldF(particles[i].getF());
            particles[i].setF(newForces[i]);
        }
    } else {
        for (auto &p1 : particles) {
            p1.setOldF(p1.getF());
            std::array<double, 3> cal_f = {0,0,0};
            for (auto &p2 : particles) {
                if(&p1 != &p2) {
                    std::array<double, 3> forceIJ = f.force(p2, p1);
                    cal_f = cal_f + forceIJ;
                }
                p1.setF(cal_f);
            }
        }
    }
}


void ParticleContainer::updateX(double delta_t){
    for (auto & particle : particles) {
        std::array<double, 3> vec = particle.getX() + delta_t * (particle.getV() + (delta_t / (2 * particle.getM())) * particle.getF());
        particle.setX(vec);
  }
}


void ParticleContainer::updateV(double delta_t){
    for (auto &p : particles) {
        std::array<double, 3> vec = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(vec);
    }
}

void ParticleContainer::simulate(double delta_t, double end_time, const std::string& out_name, const std::string& output_format, bool save_output) {

    int max_iteration = (int) (end_time / delta_t);
    for (int iteration = 0; iteration < max_iteration; iteration++) {
        if (iteration % 10 == 0 && save_output) {
            plotParticles(iteration, out_name, output_format);
        }
        // Calculate the position
        updateX(delta_t);
        // Calculate the force
        updateF();
        // Calculate the velocity
        updateV(delta_t);
        // Plot every 10th iteration
        std::cout << "Iteration " << (iteration + 1) << " finished." << "\n";
    }
    std::cout << "output written. Terminating..." << "\n";
}

void ParticleContainer::plotParticles(int iteration, const std::string& out_name, const std::string& output_format) {
    std::cout << "Plotting Particles..." << "\n";

    if(output_format == "vtu") {
        outputWriter::VTKWriter writer;
        writer.writeFile(out_name, iteration,particles);
    }
    else if(output_format == "xyz") {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particles, out_name, iteration);
    }
}

int ParticleContainer::getParticleNumber() const {
    return particles.size();
}

std::vector<Particle>& ParticleContainer::getParticles() const {
    return particles;
}

std::string ParticleContainer::toString() {
    std::stringstream buf;
    buf << "Number of particles: " << particles.size() << std::endl;
    for (auto &p : particles) {
        buf << p.toString() << std::endl;
    }
    return buf.str();
}

bool ParticleContainer::operator==(const ParticleContainer& other) const {
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
