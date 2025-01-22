cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn_GroupD VERSION 0.0.1)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.

set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)


file(GLOB_RECURSE MY_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    # header don't need to be included but this might be necessary for some IDEs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)
add_executable(MolSim ${MY_SRC})

# set c++ standard to c++20
set_target_properties(MolSim PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        ENABLE_EXPORTS ON
)

target_include_directories(MolSim
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src

)

target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
            xerces-c
            spdlog::spdlog
            fmt::fmt
)   

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -O3 -ip -ipo -fast -wd383,981,1418,1572,2259>

)


# set the build type
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp -g")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fopenmp -g")  

#set the build type to release -O3 flag is used for optimization
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp -O3 -march=native -pg")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -fopenmp -march=native -pg")

# Optionally add separate Debug and Release settings to avoid overwriting them
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in Release mode ")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in Debug mode ")
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)
include(spdlog)
include(gtest)
include (fmt)
include(xsd)
include(GoogleTest)

#testfiles 
file(GLOB_RECURSE MY_TESTS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)


list(FILTER MY_TESTS EXCLUDE REGEX "MolSim.cpp") # exclude main file from tests
add_executable(tests ${MY_TESTS}) # make tests file executable

target_include_directories(tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(tests PRIVATE GTest::gtest_main spdlog::spdlog gmock MolSim)

set_property(TARGET tests PROPERTY ENVIRONMENT "GTEST_COLOR=1")

gtest_discover_tests(tests)

enable_testing()
