/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "body/Particle.h"

#include <vector>

/**
 * @brief This class provides a function to read input data from plain .txt files.
 */
class FileReader {

public:
    FileReader();
    virtual ~FileReader();

    /**
     * Read input from plain .txt files.
     * @param particles The vector that particles should be read into.
     * @param filename The name of the input file.
     */
    void readFile(std::vector<Particle> &particles, std::string filename);
};
