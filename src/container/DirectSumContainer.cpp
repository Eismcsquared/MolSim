#include <vector>
#include <array>
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

    resetF();

    for (unsigned long i = 0; i < particles.size(); ++i) {
        for (unsigned long j = i + 1; j < particles.size(); ++j) {


            std::array<double, 3> forceIJ = force->force(particles[i], particles[j]);

            // update forces using the third Newton axiom
            particles[i].addForce(-1 * forceIJ);
            particles[j].addForce(forceIJ);
        }
    }
}


void DirectSumContainer::updateX(double delta_t){
    for (auto &particle : particles) {
        particle.updateX(delta_t);
  }
}


void DirectSumContainer::updateV(double delta_t){
    for (auto &particle : particles) {
        particle.updateV(delta_t);
    }
}





