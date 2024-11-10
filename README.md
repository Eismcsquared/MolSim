MolSim GROUP D
==============

Team Members: Siyuan Huang, Jongtae Park

## Instructions

### Build and Run

`mkdir build`<br>
`cd build`<br>
`cmake .. [-DBUILD_DOCUMENTATION=ON|OFF]`<br>
`make`<br>
`./Molsim $(pwd)/input/test.txt` #Example : `./Molsim ../input/eingabe-sonne.txt`

cmake version 3.28.3 g++ 13.2.0 doxygen 1.12.0 clang version 18.1.3

### Generate Doxygen Documentations

`make doc_doxygen`<br>

Note: Only possible if the project was built with the option
`-DBUILD_DOCUMENTATION=ON`
