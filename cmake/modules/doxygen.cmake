# make doc_doxygen optional if someone does not have / like doxygen

# TODO: create CMake build option for the target.
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Build Doxygen documentation" ON)
# TODO: Add a custom target for building the documentation.
if (BUILD_DOCUMENTATION AND ${DOXYGEN_FOUND})
    set(DOXYGEN_CONFIG "${CMAKE_SOURCE_DIR}/Doxyfile")
    add_custom_target(doc_doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} "${DOXYGEN_CONFIG}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_property(TARGET doc_doxygen PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()
