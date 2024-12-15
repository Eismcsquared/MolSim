#include <vector>
#include "body/Particle.h"

class StateReader {
public:
    static void loadState(std::vector<Particle> &particles, std::string fileName);
};

