/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "Cuboid.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(std::vector<Particle> &particles, char *filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        std::cout << "Reading " << num_particles << "." << std::endl;

        for (int i = 0; i < num_particles; i++) {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
            tmp_string.erase(tmp_string.find_last_not_of(" \t\n\r\f\v") + 1);
            std::istringstream datastream(tmp_string);
            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }
            if (datastream.eof()) {
                std::cout
                    << "Error reading file: eof reached unexpectedly reading from line "
                    << i << std::endl;
                exit(-1);
            }
            datastream >> m;
            if (datastream.eof()) {
                particles.emplace_back(x, v, m);
            } else {
                std::array<unsigned int, 3> n;
                double distance, avgV;
                for (auto &nj: n) {
                    datastream >> nj;
                }
                datastream >> distance;
                datastream >> avgV;
                Cuboid cuboid(x, v, n, m, distance, avgV);
                cuboid.createParticles(particles);
            }
        }
    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}
