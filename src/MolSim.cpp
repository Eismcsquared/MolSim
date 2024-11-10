/**
 * MolSim.cpp
 * @brief created by: Group D , main file for the simulation
 * @date 31.10.2024
 * 
 */

#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "GravitationalForce.h"

#include <iostream>
#include <list>

#include "ParticleContainer.h"
#include <cmath>
#include "outputWriter/VTKWriter.h"
#include <getopt.h>
#include <cstdlib>
#include <float.h>
#include <unistd.h>
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
    double start_time = 0;
    double end_time = 1000;
    double delta_t = 0.014;
    std::string outputFormat = ".vtu";

    std::vector<Particle> particles;
  std::cout << "Hello from MolSim for PSE!" << "\n"
  << "To run this program, please provide the file name as an argument, like this: `./Molsim abc.txt " << "\n"
  << "To see more options, type ./Molsim help or ./Molsim --help" << "\n\n";

  
  for (int i = 0; i < argc; i++) {
    if (std::string(argsv[i]) == "--help" || std::string(argsv[i]) == "help") {
      printHelp();
      return 0;
    }
  }

  if (argc < 2) {
    std::cout << "Error: Filename is required." << "\n";
    printHelp();
    return 1;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);


  if (argc >= 3) delta_t = std::atof(argsv[2]);
  if (argc >= 4) end_time = std::atof(argsv[3]);
  if (argc >= 5) {
      outputFormat = argsv[4];
    if (outputFormat != ".xyz" && outputFormat != ".vtu") {
      std::cout << "Invalid output format! Choose either 'xyz' or 'vtu'." << "\n\n";
      return 1; // exit with error
    }
  }

  // Start the timer
  auto start = std::chrono::high_resolution_clock::now();
  // Create a particle container for forwarding the particles, start time, end time, delta_t and outputFormat
  GravitationalForce f;
  ParticleContainer particle_container = ParticleContainer(particles, start_time, end_time, delta_t, f, outputFormat);

  // Inform the user about the input parameters
  std::cout << "Testfilename: " << argsv[1] << "\n";
  std::cout << "Start Time: " << start_time << "\n";
  std::cout << "Time End: " << end_time << "\n" ;
  std::cout << "Delta Time: " << delta_t << "\n";
  std::cout << "Output format: " << outputFormat << "\n\n";


  // Calculate the position, force and velocity for all particles
    particle_container.simulate();

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