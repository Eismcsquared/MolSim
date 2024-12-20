/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "body/Cuboid.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "spdlog/spdlog.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(std::vector<Particle> &particles, std::string filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);

        spdlog::trace("Read line: {}", tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            spdlog::trace("Read line: {}", tmp_string);
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        spdlog::trace("Reading {}.", num_particles);

        for (int i = 0; i < num_particles; i++) {
            getline(input_file, tmp_string);
            spdlog::trace("Read line: {}", tmp_string);
            tmp_string.erase(tmp_string.find_last_not_of(" \t\n\r\f\v") + 1);
            std::istringstream datastream(tmp_string);
            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }
            if (datastream.eof()) {
                spdlog::error("Error reading file: eof reached unexpectedly reading from line {}", i);
                std::exit(-1);
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
                Cuboid cuboid(x, v, n, m, distance, avgV, 3);
                cuboid.createParticles(particles);
            }
        }
    } else {
        spdlog::error("Error: could not open file {}", filename);
        exit(-1);
    }
}
