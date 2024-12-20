#include <fstream>
#include <vector>
#include "body/Particle.h"

/**
 * @brief This class saves the current state of a system to a file for checkpointing.
 */
class StateWriter {
public:

    static void saveState(std::vector<Particle> &particles, std::string fileName);
};
