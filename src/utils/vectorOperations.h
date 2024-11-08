#include <array>
#include <cmath>

std::array<double, 3> add(std::array<double, 3> array1, std::array<double, 3> array2) {
    return {array1[0] + array2[0], array1[1] + array2[1], array1[2] + array2[2]};
}

std::array<double, 3> sub(std::array<double, 3> array1, std::array<double, 3> array2) {
    return {array1[0] - array2[0], array1[1] - array2[1], array1[2] - array2[2]};
}

std::array<double, 3> scalar_mul(double c, std::array<double, 3> array) {
    return {c * array[0], c * array[1], c * array[2]};
}

double norm(std::array<double, 3> array) {
    return sqrt(pow(array[0], 2) + pow(array[1], 2) + pow(array[2], 2));
}
