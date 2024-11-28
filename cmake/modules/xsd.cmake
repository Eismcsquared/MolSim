find_program(XSD_EXECUTABLE xsd)

if (XSD_EXECUTABLE)
    file(GLOB_RECURSE XSD
            "${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader/*.xsd"
    )
    add_custom_target(xsd
            COMMAND ${XSD_EXECUTABLE} cxx-tree --cxx-flags "-std=c++17" --hxx-suffix .h --cxx-suffix .cpp --generate-doxygen --generate-serialization ${XSD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader
            COMMENT "Generating code from xml schemas"
            VERBATIM
    )

endif()