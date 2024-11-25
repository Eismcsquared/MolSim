#include <algorithm>
#include "Cell.h"

Cell::Cell(std::array<double, 3> position, std::array<double, 3> size): position(position), size(size) {}

const std::vector<Particle> Cell::getParticles() {
    return particles;
}

const std::array<double, 3> Cell::getPosition() {
    return position;
}

const std::array<double, 3> Cell::getSize() {
    return size;
}

bool Cell::contains(std::array<double, 3> pos) {
    for (int i = 0; i < 3; ++i) {
        if (pos[i] < position[i] || pos[i] > position[i] + size[i]) {
            return false;
        }
    }
    return true;
}

void Cell::addParticle(Particle &particle) {
    if (contains(particle.getX())) {
        particles.insert(particles.end(), particle);
    }
}

void Cell::removeParticle(Particle &particle) {
    auto newEnd = std::remove(particles.begin(), particles.end(), particle);
    particles.erase(newEnd, particles.end());
}
