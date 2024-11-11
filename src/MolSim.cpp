/**
 * MolSim.cpp
 * @brief created by: Group D , main file for the simulation
 * @date 31.10.2024
 * 
 */

#include "FileReader.h"
#include "GravitationalForce.h"
#include "LennardJonesForce.h"

#include <iostream>
#include <list>

#include "ParticleContainer.h"
#include <getopt.h>
#include <cstdlib>
#include <vector>
#include <chrono>


/**
 * @brief Function to plot the particles
 */
void plotParticles(int iteration);

/**
 * @brief Function to print the help message
 */
void printHelp();



int main(int argc, char *argsv[]) {
    enum ForceType {
        GRAVITATION,
        LENNARD_JONES
    };
    double start_time = 0;
    double end_time = 5;
    double delta_t = 0.0002;
    std::string outputFormat = "vtu";
    ForceType mode = GRAVITATION;
    if (argc < 2) {
        std::cout << "Error: Filename is required." << "\n";
        printHelp();
        return 1;
    }
    char* inputFile = argsv[1];
    std::string outputFile("MD_vtk");

    std::vector<Particle> particles;
    std::cout << "Hello from MolSim for PSE!" << "\n"
    << "To run this program, please provide the input file name as an argument, like this: `./Molsim abc.txt " << "\n"
    << "To see more options, type ./Molsim help or ./Molsim --help" << "\n\n";


    int opt;
    static struct option long_options[] = {
            {"help",    no_argument,       nullptr,  'h' },
            {"format", required_argument, nullptr, 'f'},
            {"end_time", required_argument, nullptr, 'e'},
            {"delta_t", required_argument, nullptr, 'd'},
            {"output", required_argument, nullptr, 'o'},
            {"gravitation", no_argument, nullptr, 'g'},
            {"Lennard_Jones", no_argument, nullptr, 'l'}

    };
    while ((opt = getopt_long(argc, argsv, "hd:e:f:o:gl", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'd':
                delta_t = std::atof(optarg);
                break;
            case 'e':
                end_time = std::atof(optarg);
                break;
            case 'f':
                outputFormat = optarg;
                if (outputFormat != "xyz" && outputFormat != "vtu") {
                    std::cout << "Invalid output format! Choose either 'xyz' or 'vtu'." << "\n\n";
                    return 1; // exit with error
                }
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'g':
                mode = GRAVITATION;
                break;
            case 'l':
                mode = LENNARD_JONES;
                break;
            case '?':
                std::cout << "Invalid option!" << "\n";
                printHelp();
                return 1;
            case 'h':
            default:
                printHelp();
                return 0;
        }
    }

    for(int i = 0; i < argc; i++) {
        std::cout << argsv[i] << "\n";
    }

    FileReader fileReader;
    fileReader.readFile(particles, inputFile);

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();
    // Create a particle container for forwarding the particles, start time, end time, delta_t and outputFormat
    Force* force;
    switch (mode) {
        case GRAVITATION:
            force = new GravitationalForce();
            break;
        case LENNARD_JONES:
            force = new LennardJonesForce();
    }
    ParticleContainer particle_container = ParticleContainer(particles, start_time, end_time, delta_t, force, outputFormat);

  // Inform the user about the input parameters
    std::cout << "Testfilename: " << inputFile << "\n";
    std::cout << "Start Time: " << start_time << "\n";
    std::cout << "Time End: " << end_time << "\n" ;
    std::cout << "Delta Time: " << delta_t << "\n";
    std::cout << "Output format: " << outputFormat << "\n\n";


    // Calculate the position, force and velocity for all particles
    particle_container.simulate(outputFile);

    std::cout << "output written. Terminating..." << "\n";
   
  
    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

  
    // Calculate the duration
    std::chrono::duration<double> duration = end - start;
    std::cout << "Duration: " << duration.count() << "s" << "\n";


    return 0;
}


void printHelp() {
    std::cout << "MolSim for PSE - Usage:" << "\n";
    std::cout << "./MolSim filename [delta_t] [end_time] [output_format]" << "\n";
    std::cout << "Options:" << "\n";
    std::cout << "  1.  Specify the input file name (required)" << "\n";
    std::cout << "  2.  Set the delta time (optional, default is 0.014)" << "\n";
    std::cout << "  3.  Set the end time (optional, default is 1000)" << "\n";
    std::cout << R"(  4.  Specify output format: either ".xyz" or ".vtu" (optional, default is ".vtu "))" << "\n";
    std::cout << "  --help, help message" << "\n\n";
}