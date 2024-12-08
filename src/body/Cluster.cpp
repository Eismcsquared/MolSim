#include "Cluster.h"

Cluster::Cluster(const std::array<double, 3> &x, const std::array<double, 3> &v, double m,
                 double distance, double avgVelocityBrownian, int dimension):
                 x(x), v(v), m(m), distance(distance), avgVelocityBrownian(avgVelocityBrownian), dimension(dimension),
                 epsilon(5), sigma(1){}

Cluster::Cluster(const std::array<double, 3> &x, const std::array<double, 3> &v, double m,
                 double distance, double avgVelocityBrownian, int dimension, double epsilon,
                 double sigma):
        x(x), v(v), m(m), distance(distance), avgVelocityBrownian(avgVelocityBrownian), dimension(dimension),
        epsilon(epsilon), sigma(sigma){}
