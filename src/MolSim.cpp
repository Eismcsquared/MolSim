/**
 * MolSim.cpp
 * @brief created by: Group D , main file for the simulation
 * @date 31.10.2024
 * 
 */

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <memory>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "inputReader/FileReader.h"
#include "inputReader/XMLReader.h"
#include "force/GravitationalForce.h"
#include "force/LennardJonesForce.h"
#include "container/DirectSumContainer.h"
#include "container/LinkedCellContainer.h"
#include "Simulation.h"
/**
 * @brief Function to plot the particles
 */
void plotParticles(int iteration);

/**
 * @brief Function to print the help message
 */
void printHelp();


int main(int argc, char *argsv[]) {
    spdlog::set_level(spdlog::level::info);
    spdlog::stdout_color_mt("console"); // Create color multithreaded logger

    spdlog::trace("MolSim started");

    if (argc < 2) {
        spdlog::error("Error: Filename is required.");
        printHelp();
        return 1;
    }

    std::string inputFile= std::string (argsv[1]);
    std::unique_ptr<Simulation> simulation;

    if (inputFile.size() >= 4 && inputFile.compare(inputFile.size() - 4, 4, ".xml") == 0) {
        simulation = XMLReader::readXML(inputFile);
    } else {
        // default values when using .txt as input
        double start_time = 0;
        double end_time = 5;
        double delta_t = 0.0002;
        unsigned int frequency = 10;
        std::string outputFile("MD_vtk");
        std::string outputFormat = "vtu";

        std::unique_ptr<std::vector<Particle>> particles = std::make_unique<std::vector<Particle>>();
        FileReader fileReader;
        fileReader.readFile(*particles, inputFile);
        std::unique_ptr<Force> f = std::make_unique<LennardJonesForce>();
        std::unique_ptr<ParticleContainer> container = std::make_unique<DirectSumContainer>(particles, f);
        simulation = std::make_unique<Simulation>(container, end_time, delta_t, outputFile, outputFormat, frequency);
    }

    std::unique_ptr<Force> f;
    int opt;
    static struct option long_options[] = {
            {"help",    no_argument,       nullptr,  'h' },
            {"format", required_argument, nullptr, 'f'},
            {"end_time", required_argument, nullptr, 'e'},
            {"delta_t", required_argument, nullptr, 'd'},
            {"output", required_argument, nullptr, 'o'},
            {"spdlog_level", required_argument, nullptr, 's'},
            {"benchmark", no_argument, nullptr, 'b'},
            {"newton3 ", no_argument, nullptr, 'n'},
            {"gravitation", no_argument, nullptr, 'g'},
            {"Lennard_Jones", no_argument, nullptr, 'l'},
            {nullptr, 0, nullptr, 0} 
    };
    while ((opt = getopt_long(argc, argsv, "hd:e:f:o:s:glbn", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'd':
                simulation->setDeltaT(std::atof(optarg));
                break;
            case 'e':
                simulation->setEndTime(std::atof(optarg));
                break;
            case 'f':
                if ((std::string )optarg != "xyz" && (std::string )optarg != "vtu") {
                    spdlog::error("Invalid output format! Choose either 'xyz' or 'vtu'.");
                    return 1; // exit with error
                }
                simulation->setOutputFormat((std::string) optarg);
                break;
            case 'o':
                simulation->setOutputFile((std::string) optarg);
                break;
            case 's':
                if (optarg) {
                    int lv = std::atoi(optarg);  
                    switch (lv) {
                        case 1: spdlog::set_level(spdlog::level::trace); break;
                        case 2: spdlog::set_level(spdlog::level::debug);  break;
                        case 3: spdlog::set_level(spdlog::level::info); break;
                        case 4: spdlog::set_level(spdlog::level::warn); break;
                        case 5: spdlog::set_level(spdlog::level::err); break;
                        case 6: spdlog::set_level(spdlog::level::trace); break;
                        default:
                            spdlog::error("Invalid spdlog level! Choose a number from 1 to 6.");
                            return 1;
                    }
                } else {  
                    if(optarg == nullptr){
                        spdlog::error("No log levasel specified for -s option.");
                        return 1;
                    }
                    spdlog::error("No log level specified for -s option.");
                    return 1;
                }
                break;
            case 'b':
                simulation->setSaveOutput(false);
                spdlog::set_level(spdlog::level::off);
                break;
            case 'n':
                simulation->setNewton3(false);
                break;
            case 'g':
                f = std::make_unique<GravitationalForce>();
                simulation->getContainer()->setF(f);
                break;
            case 'l':
                f = std::make_unique<LennardJonesForce>();
                simulation->getContainer()->setF(f);
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

  // Inform the user about the input parameters
    spdlog::info("Input file name: {}", inputFile);
    spdlog::info("Start Time: {}", 0);
    spdlog::info("Time End: {}", simulation->getEndTime());
    spdlog::info("Delta Time: {}", simulation->getDeltaT());
    spdlog::info("Output file name: {}", simulation->getOutputFile());
    spdlog::info("Output format: {}", simulation->getOutputFormat());

    auto start = std::chrono::high_resolution_clock::now();
    // run simulation
    simulation->run();
    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    std::chrono::duration<double> duration = end - start;


    spdlog::info("output written. Terminating...");
    spdlog::set_level(spdlog::level::info);

    spdlog::info("Duration: {}s", duration.count());

    return 0;
}


void printHelp() {
    std::cout << "MolSim for PSE" << "\n";
    std::cout << "Usage:" << "\n";
    std::cout << "  ./MolSim <input-file> [-d <time-step>] [-e <duration>] [-f <output-format>] [-o <output-file>] [-l|-g] [-s <log-level>]" << "\n";
    std::cout << "Options:" << "\n";
    std::cout << "  -d or --delta_t <time-step>    = The length of each time step of the simulation (Default: 0.0002)." << "\n";
    std::cout << "  -e or --end_time <duration>    = The total duration of the simulation (Default = 5)." << "\n";
    std::cout << "  -f or --format <output-format> = The format of the output, must be either vtu or xyz (Default: vtu)." << "\n";
    std::cout << "  -o or --output <output-file>   = The name of files that data should be written to (Default: MD_vtk)." << "\n";
    std::cout << "  -s or --spdlog_level <level>   = Set spdlog level (trace -1, debug -2, info -3, warn -4, error -5, critical -6).\n";
    std::cout << "  -b or --benchmark              = If specified, the benchmark mode is activated." << "\n";
    std::cout << "  -n or --newton3                = If specified, the Newton's third law is not applied." << "\n";
    std::cout << "  -g or --gravitation            = The simulation of gravitational force (with G = 1) between objects (Default)." << "\n";
    std::cout << "  -l or --Lennard_Jones          = If specified, the Lennard Jones potential (with epsilon = 5 and sigma = 1) is simulated." << "\n";
    std::cout << "  -h or --help                   = Print help message." << "\n";
}