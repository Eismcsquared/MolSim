#include "Sphere.h"


Sphere::Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius, double m,
               double distance, double avgVelocityBrownian, int dimension) :
               Cluster(x, v, m, distance, avgVelocityBrownian, dimension), radius(radius){}


Sphere::Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius, double m, double distance,
               double avgVelocityBrownian, int dimension, double epsilon, double sigma):
        Cluster(x, v, m, distance, avgVelocityBrownian, dimension, epsilon, sigma), radius(radius){

}


