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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

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

    auto logger = spdlog::stdout_color_mt("console");
    logger->set_level(spdlog::level::trace);
    spdlog::trace("MolSim started");
    char* filename = const_cast<char*>("../input/assignment2.txt");


    if (argc < 2) {
        spdlog::error("Error: Filename is required.");
        printHelp();
        return 1;
    }

    char* inputFile = filename;
    inputFile= argsv[1];
    std::string outputFile("MD_vtk");

    std::vector<Particle> particles;

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
                    spdlog::error("Invalid output format! Choose either 'xyz' or 'vtu'.");
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
                spdlog::error("Invalid option!");
                printHelp();
                return 1;
            case 'h':
            default:
                printHelp();
                return 0;
        }
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
    ParticleContainer particle_container = ParticleContainer(particles, force);

  // Inform the user about the input parameters
    spdlog::info("Testfilename: {}", inputFile);
    spdlog::info("Start Time: {}", start_time);
    spdlog::info("Time End: {}", end_time);
    spdlog::info("Delta Time: {}", delta_t);
    spdlog::info("Output format: {}", outputFormat);
    

    // Calculate the position, force and velocity for all particles
    particle_container.simulate(delta_t, end_time, outputFile, outputFormat);

    spdlog::info("output written. Terminating...");
   
  
    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

  
    // Calculate the duration
    std::chrono::duration<double> duration = end - start;
    spdlog::info("Duration: {}s", duration.count());


    return 0;
}


void printHelp() {
    std::cout << "MolSim for PSE" << "\n";
    std::cout << "Usage:" << "\n";
    std::cout << "  ./MolSim <input-file> [-d <time-step>] [-e <duration>] [-f <output-format>] [-o <output-file>] [-l|-g]" << "\n";
    std::cout << "Options:" << "\n";
    std::cout << "  -d or --delta_t <time-step>    = The length of each time step of the simulation (Default: 0.0002)." << "\n";
    std::cout << "  -e or --end_time <duration>    = The total duration of the simulation (Default = 5)." << "\n";
    std::cout << "  -f or --format <output-format> = The format of the output, must be either vtu or xyz (Default: vtu)." << "\n";
    std::cout << "  -o or --output <output-file>   = The name of files that data should be written to (Default: MD_vtk)." << "\n";
    std::cout << "  -g or --gravitation            = The simulation of gravitational force (with G = 1) between objects (Default)." << "\n";
    std::cout << "  -l or --Lennard_Jones          = If specified, the Lennard Jones potential (with epsilon = 5 and sigma = 1) is simulated." << "\n";
    std::cout << "  -h or --help                   = Print help message." << "\n";
}