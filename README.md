MolSim GROUP D
==============

Team Members: Siyuan Huang, Jongtae Park

## Instructions

### Requirements

cmake version 3.28.3, g++ 13.2.0, doxygen 1.12.0, clang version 18.1.3

### Build and Run

`mkdir build`<br>
`cd build`<br>
`cmake .. [-DBUILD_DOCUMENTATION=ON|OFF]`<br>
`make`<br>
`./Molsim <input-file>` #Example : `./Molsim ../input/eingabe-sonne.txt`

### Usage

```
./MolSim <input-file> [-d <time-step>] [-e <duration>] [-f <output-format>] [-o <output-file>] [-l|-g]
```

Options:

`-d | --delta_t <time-step>` = The length of each time step of the simulation (Default: 0.0002).<br>
`-e | --end_time <duration>` = The total duration of the simulation (Default = 5).<br>
`-f | --format <output-format>` = The format of the output, must be either vtu or xyz (Default: vtu).<br>
`-o | --output <output-file>` = The name of files that data should be written to (Default: MD_vtk).<br>
`-s | --spdlog_level <level>` = Set spdlog level (trace -1, debug -2, info -3, warn -4, error -5, critical -6).<br>
`-g | --gravitation` = The simulation of gravitational force (with G = 1) between objects.<br>
`-l | --Lennard_Jones` = If specified, the Lennard Jones potential (with epsilon = 5 and sigma = 1) is simulated.<br>
`-h | --help` = Print help message.<br>


Input file format:<br>

The first line specify the number of data sets. Each of the following line is interpreted as a data set, that either represents a single particle or a cuboid of particles. Each data set consists of the following data, seperated by white spaces:<br>

- The spatial coordinates (x, y, z) of the object, seperated by white spaces
- The velocity (v<sub>x</sub>, v<sub>y</sub>, v<sub>z</sub>) of the objects, seperated by white spaces
- The mass of the/each particle

For a cuboid of particles, The spatial coordinates specified above are interpreted as the position of the particle with the least coordinates and the velocity is the average velocity of particles in the cuboid. In addition, the following parameters should be specified:

- The number of particles in each direction (N<sub>x</sub>, N<sub>y</sub>, N<sub>z</sub>), separated by white spaces
- The distance between two neighbouring particles
- The average velocity of the Brownian motion

Comments are mark with "#" at the beginning of a line and are only allowed at the beginning of the file.

### Generate Doxygen Documentations

`make doc_doxygen`<br>

Note: Only possible if the project was built with the option
`-DBUILD_DOCUMENTATION=ON`

### Assignments

For simulations required in the individuell work sheets, run the following commands.<br>

- Assignment 1:
  `./MolSim ../input/eingabe-sonne.txt -d 0.014 -e 1000 -g`<br>
- Assignment 2:
  `./MolSim ../input/assignment2.txt -d 0.0002 -e 5 -l`<br>

