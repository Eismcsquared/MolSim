#include <algorithm>
#include "Cell.h"

Cell::Cell(std::array<double, 3> position, std::array<double, 3> size, std::set<int>& neighbours):
position(position), size(size), neighbours(neighbours) {}

const std::vector<int> &Cell::getParticleIndices() const {
    return particleIndices;
}

const std::array<double, 3> Cell::getPosition() const {
    return position;
}

const std::array<double, 3> Cell::getSize() const {
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
    particleIndices.insert(particleIndices.end(), index);
}

void Cell::removeIndex(int index) {
    auto newEnd = std::remove(particleIndices.begin(), particleIndices.end(), index);
    particleIndices.erase(newEnd, particleIndices.end());
}

void Cell::clear() {
    particleIndices.clear();
}
