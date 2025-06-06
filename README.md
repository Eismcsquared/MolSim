MolSim GROUP D
==============

Team Members: Siyuan Huang, Jongtae Park

## Instructions

### Requirements

cmake version 3.28.3, g++ 13.2.0, doxygen 1.12.0, clang version 18.1.3

### Build and Run

`mkdir build`<br>
`cd build`<br>
`cmake .. [-DBUILD_DOCUMENTATION=ON|OFF] [-DOPENMP=ON|OFF]`<br>
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
`-b | --benchmark` = If specified, the benchmark mode is activated and the output of data is deactivated.<br>
`-g | --gravitation` = If specified, The gravitation (with G = 1) is taken to be the force between objects. Otherwise, the Lennard-Jones force is assumed by default.<br>
`-c | --checkpoint <file>` = If specified, the final state is stored to the given file and the output of simulation data is deactivated.<br>
`-l | --load <file>` = If specified, additional particles are loaded from the given file.<br>
`-h | --help` = Print help message.<br>


### Input file format

#### Plain text input

The first line specify the number of data sets. Each of the following line is interpreted as a data set, that either represents a single particle or a cuboid of particles. Each data set consists of the following data, seperated by white spaces:<br>

- The spatial coordinates (x, y, z) of the object, seperated by white spaces
- The velocity (v<sub>x</sub>, v<sub>y</sub>, v<sub>z</sub>) of the objects, seperated by white spaces
- The mass of the/each particle

For a cuboid of particles, The spatial coordinates specified above are interpreted as the position of the particle with the least coordinates and the velocity is the average velocity of particles in the cuboid. In addition, the following parameters should be specified:

- The number of particles in each direction (N<sub>x</sub>, N<sub>y</sub>, N<sub>z</sub>), separated by white spaces
- The distance between two neighbouring particles
- The average velocity of the Brownian motion

Comments are mark with "#" at the beginning of a line and are only allowed at the beginning of the file.

#### XML file input

The XML input scheme is defined in `src/inputReader/InputData.xsd`. Note that the plain text input format does not support all parameters that are supported by the XML scheme. Hence, the XML input format is preferred.

### Generate Doxygen Documentations

`make doc_doxygen`<br>

Note: Only possible if the project was built with the option
`-DBUILD_DOCUMENTATION=ON`

### Run tests
`ctest` or `./tests`

### Parallelization
`cmake -DOPENMP=ON ..`<br>
`make`<br>
`[OMP_NUM_THREADS=<number-of-threads>] ./MolSim <input-file>`

### Assignments

For simulations required in the individuell work sheets, run the following commands.<br>

- Assignment 1:
  `./MolSim ../input/assignment1.xml`<br>
- Assignment 2:
  `./MolSim ../input/assignment2.xml`<br>
- Assignment 3 - Collision:
  `./MolSim ../input/assignment3.xml`
- Assignment 3 - Falling drop:
  `./MolSim ../input/assignment3_falling_drop.xml`
- Assignment 4 - Rayleigh-Taylor instability: 
  `./MolSim ../input/assignment4.xml`
- Assignment 4 - Falling drop:
  `./MolSim ../input/assignment4_fluid.xml` followed by
  `./MolSim ../input/assignment4_falling_drop.xml`
- Assignment 5 - Membrane:
  `./MolSim ../input/assignment5_membrane.xml`
- Assignment 5 - Rayleigh-Taylor instability:
  `./MolSim ../input/assignment5.xml`
- Assignment 5 - Nano-scale Flow:
  `./MolSim ../input/assignment5_flow.xml`

### Benchmark
Comparison: Linked cell algorithm vs. Direct sum algorithm for different number of particles. 
![](images/ds_vs_lc.png)

Comparison: Run time and molecular updates per second before and after optimizations. The measurements are based on the input assignment4.xml and run on Linux cluster.
![](images/assignment4_task2.png)

Comparison: Speedup for Varying Thread Numbers. The strong scaling test was performed for the thread count 1, 2, 4, 8, 14, 16, 28, and 56 over 1000 iterations. The measurements were based on the input file assignment5.xml and run on the Linux cluster cm4_tiny. 
![](images/Speadup_gcc.png)
