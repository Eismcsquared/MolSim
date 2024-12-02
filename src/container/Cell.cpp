#include <algorithm>
#include "Cell.h"

Cell::Cell(std::array<double, 3> position, std::array<double, 3> size): position(position), size(size) {}

const std::vector<unsigned int> &Cell::getParticleIndices() {
    return particleIndices;
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

void Cell::addIndex(unsigned int index) {
    particleIndices.insert(particleIndices.end(), index);
}

void Cell::removeIndex(unsigned int index) {
    auto newEnd = std::remove(particleIndices.begin(), particleIndices.end(), index);
    particleIndices.erase(newEnd, particleIndices.end());
}
