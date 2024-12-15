#include <iostream>
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
    // compute initial forces
    ParticleContainer::updateF();
    spdlog::trace("DirectSumContainer generated!");
}


DirectSumContainer::~DirectSumContainer(){
    spdlog::trace("DirectSumContainer destructed!");
}


void DirectSumContainer::updateF(bool newton3) {
    if (newton3) {
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
    } else {
        for (auto &p1 : particles) {
            p1.setOldF(p1.getF());
            std::array<double, 3> cal_f = {0, p1.getM() * g, 0};
            for (auto &p2 : particles) {
                if(&p1 != &p2) {
                    std::array<double, 3> forceIJ = force->force(p2, p1);
                    cal_f = cal_f + forceIJ;
                }

            }
            p1.setF(cal_f);
        }
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





