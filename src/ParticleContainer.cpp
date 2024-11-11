#include <iostream>
#include <vector>
#include <array>
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include "ParticleContainer.h"

/**
 * @brief Constructor, which initializes the particle container with a vector of particles , start time, end time, delta_t and outputFormat
 * outputFormat can be either .vtu or .xyz
 */
ParticleContainer::ParticleContainer(std::vector<Particle>& particles, double start_time,
                                     double end_time, double delta_t, Force& f, std::string outputFormat)
    : particles(particles),
      start_time(start_time),
      end_time(end_time),
      delta_t(delta_t),
      f(f),
      outputFormat(std::move(outputFormat)) { // std::move if outputFormat is temporary
    // compute initial forces
    updateF();
    std::cout << "ParticleContainer generated!\n";
}



ParticleContainer::~ParticleContainer(){
    std::cout << "ParticleContainer destructed!\n";
}


void ParticleContainer::addParticle(const Particle& particle){
    this->particles.push_back(particle);
}

void ParticleContainer::addCuboid(const Cuboid& cuboid) {
    cuboid.createParticles(particles);
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


void ParticleContainer::updateX(){
    for (auto & particle : particles) {
        std::array<double, 3> vec = particle.getX() + delta_t * (particle.getV() + (delta_t / (2 * particle.getM())) * particle.getF());
        particle.setX(vec);
  }
}


void ParticleContainer::updateV(){
    for (auto &p : particles) {
        std::array<double, 3> vec = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(vec);
    }
}

void ParticleContainer::simulate() {

    int max_iteration = (int) ((end_time - start_time) / delta_t);
    for (int iteration  = 0; iteration <= max_iteration; iteration++) {
        // Calculate the position
        updateX();
        // Calculate the force
        updateF();
        // Calculate the velocity
        updateV();
        // Plot every 10th iteration
        if (iteration % 10 == 0) {
            plotParticles(iteration);
        }
        std::cout << "Iteration " << iteration << " finished." << "\n";
    }
    std::cout << "output written. Terminating..." << "\n";
}

void ParticleContainer::plotParticles(int iteration) {
    std::cout << "Plotting Particles..." << "\n";

    std::string out_name("MD_vtk");
    if(outputFormat.compare(".vtu") == 0) {
        outputWriter::VTKWriter writer2;
        writer2.writeFile(out_name, iteration,particles);
    }
    else if(outputFormat.compare(".xyz") == 0) {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particles, out_name, iteration);
    }
}
