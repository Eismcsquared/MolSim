#include "Cuboid.h"

Cuboid::Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
               double m, double distance, double avgVelocityBrownian, int dimension) :
               Cluster(x, v, m, distance, avgVelocityBrownian, dimension), N(n) {}

Cuboid::Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
               double m, double distance, double avgVelocityBrownian, int dimension, int type, double epsilon, double sigma):
        Cluster(x, v, m, distance, avgVelocityBrownian, dimension, type, epsilon, sigma), N(n){

}
