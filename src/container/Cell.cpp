#include <algorithm>
#include "Cell.h"

Cell::Cell(std::array<double, 3> position, std::array<double, 3> size, std::set<int>& neighbours):
position(position), size(size), neighbours(neighbours) {
#ifdef _OPENMP
    omp_init_lock(&monitor);
#endif
}

Cell::Cell(const Cell &other):
        position(other.position), size(other.size), particleIndices(other.particleIndices), neighbours(other.neighbours) {
#ifdef _OPENMP
    omp_init_lock(&monitor);
#endif
}

Cell::~Cell() {
#ifdef _OPENMP
    omp_destroy_lock(&monitor);
#endif
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
#ifdef _OPENMP
    omp_set_lock(&monitor);
    particleIndices.insert(particleIndices.end(), index);
    omp_unset_lock(&monitor);
#else
    particleIndices.insert(particleIndices.end(), index);
#endif
}

void Cell::removeIndex(int index) {
#ifdef _OPENMP
    omp_set_lock(&monitor);
    auto newEnd = std::remove(particleIndices.begin(), particleIndices.end(), index);
    particleIndices.erase(newEnd, particleIndices.end());
    omp_unset_lock(&monitor);
#else
    auto newEnd = std::remove(particleIndices.begin(), particleIndices.end(), index);
    particleIndices.erase(newEnd, particleIndices.end());
#endif
}

void Cell::clear() {
#ifdef _OPENMP
    omp_set_lock(&monitor);
    particleIndices.clear();
    omp_unset_lock(&monitor);
#else
    particleIndices.clear();
#endif
}
