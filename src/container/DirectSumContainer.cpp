#include <vector>
#include <array>
#include <filesystem>
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "body/Particle.h"
#include "DirectSumContainer.h"
#include "spdlog/spdlog.h"


DirectSumContainer::DirectSumContainer(std::vector<Particle>& particles, std::unique_ptr<Force> &f) : ParticleContainer(particles, f)
{

    spdlog::trace("DirectSumContainer generated!");
}


DirectSumContainer::~DirectSumContainer(){
    spdlog::trace("DirectSumContainer destructed!");
}


void DirectSumContainer::updateF() {

    std::vector<std::array<double, 3>> newForces;
    newForces.reserve(particles.size());
    for (const auto & particle : particles) {
        newForces.push_back({0, particle.getM() * g, 0});
    }
    for (unsigned long i = 0; i < particles.size(); ++i) {
        for (unsigned long j = i + 1; j < particles.size(); ++j) {


            std::array<double, 3> forceIJ = force->force(particles[i], particles[j]);

            // update forces using the third Newton axiom
            newForces[i] = newForces[i] - forceIJ;
            newForces[j] = newForces[j] + forceIJ;
        }
    }
    for (unsigned long i = 0; i < particles.size(); ++i) {
        particles[i].setOldF(particles[i].getF());
        particles[i].setF(newForces[i]);
    }
}


void DirectSumContainer::updateX(double delta_t){
    for (auto & particle : particles) {
        std::array<double, 3> vec = particle.getX() + delta_t * (particle.getV() + (delta_t / (2 * particle.getM())) * particle.getF());
        particle.setX(vec);
  }
}


void DirectSumContainer::updateV(double delta_t){
    for (auto &p : particles) {
        std::array<double, 3> vec = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(vec);
    }
}





