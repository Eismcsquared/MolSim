#include "Membrane.h"

Membrane::Membrane(std::array<double, 3> &x, std::array<double, 3> &v, std::array<unsigned int, 2> &n, double m,
                   double distance, double avgVelocityBrownian, int dimension, double k, double r0):
        Cuboid(x, v, {n[0], n[1], 0}, m, distance, avgVelocityBrownian, dimension), k(k), r0(r0){

}

void Membrane::createParticles(std::vector<Particle> &particles) const {
    unsigned long oldSize = particles.size();
    Cuboid::createParticles(particles);
    for (unsigned long i = oldSize; i < particles.size(); ++i) {
        particles[i].setK(k);
        particles[i].setR0(r0);
        int idxX = (i - oldSize) / N[1];
        int idxY = (i - oldSize) % N[1];
        if (idxX > 0) {
            particles[i].addNeighbour(particles[i - N[1]]);
            if (idxY > 0) {
                particles[i].addDiagonalNeighbour(particles[i - N[1] - 1]);
            }
            if (idxY < N[1] - 1) {
                particles[i].addDiagonalNeighbour(particles[i - N[1] + 1]);
            }
        }
        if (idxX < N[0] - 1) {
            particles[i].addNeighbour(particles[i + N[1]]);
            if (idxY > 0) {
                particles[i].addDiagonalNeighbour(particles[i + N[1] - 1]);
            }
            if (idxY < N[1] - 1) {
                particles[i].addDiagonalNeighbour(particles[i + N[1] + 1]);
            }
        }
        if (idxY > 0) {
            particles[i].addNeighbour(particles[i - 1]);
        }
        if (idxY < N[1] - 1) {
            particles[i].addNeighbour(particles[i + 1]);
        }
    }
}


