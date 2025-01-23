#include <algorithm>
#include "Cell.h"

Cell::Cell(std::array<double, 3> position, std::array<double, 3> size, std::set<int>& neighbours):
position(position), size(size), neighbours(neighbours) {
    omp_init_lock(&monitor);
}

Cell::Cell(const Cell &other):
        position(other.position), size(other.size), particleIndices(other.particleIndices), neighbours(other.neighbours) {
    omp_init_lock(&monitor);
}

Cell::~Cell() {
    omp_destroy_lock(&monitor);
}

const std::vector<int> &Cell::getParticleIndices() const {
    return particleIndices;
}

const std::array<double, 3> &Cell::getPosition() const {
    return position;
}

const std::array<double, 3> &Cell::getSize() const {
    return size;
}

const std::set<int> &Cell::getNeighbours() const {
    return neighbours;
}

bool Cell::contains(std::array<double, 3> pos) {
    for (int i = 0; i < 3; ++i) {
        if (pos[i] < position[i] || pos[i] > position[i] + size[i]) {
            return false;
        }
    }
    return true;
}

void Cell::addIndex(int index) {
    omp_set_lock(&monitor);
    particleIndices.insert(particleIndices.end(), index);
    omp_unset_lock(&monitor);
}

void Cell::removeIndex(int index) {
    omp_set_lock(&monitor);
    auto newEnd = std::remove(particleIndices.begin(), particleIndices.end(), index);
    particleIndices.erase(newEnd, particleIndices.end());
    omp_unset_lock(&monitor);
}

void Cell::clear() {
    particleIndices.clear();
}
