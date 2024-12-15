#include "Cluster.h"

Cluster::Cluster(const std::array<double, 3> &x, const std::array<double, 3> &v, double m,
                 double distance, double avgVelocityBrownian, int dimension):
                 x(x), v(v), m(m), distance(distance), avgVelocityBrownian(avgVelocityBrownian), dimension(dimension), type(0),
                 epsilon(5), sigma(1){}

Cluster::Cluster(const std::array<double, 3> &x, const std::array<double, 3> &v, double m,
                 double distance, double avgVelocityBrownian, int dimension, int type, double epsilon,
                 double sigma):
        x(x), v(v), m(m), distance(distance), avgVelocityBrownian(avgVelocityBrownian), dimension(dimension), type(type),
        epsilon(epsilon), sigma(sigma){}
