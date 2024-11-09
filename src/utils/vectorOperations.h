#include <array>
#include <cmath>


double distance(std::array<double, 3> p1, std::array<double, 3> p2) {
    return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}