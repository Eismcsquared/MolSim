#include <vector>
#include "body/Particle.h"

/**
 * @brief The class loads the state of a system that is stored by StateWriter.
 */
class StateReader {
public:
    static void loadState(std::vector<Particle> &particles, std::string fileName);
};

