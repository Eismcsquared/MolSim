find_program(XSD_EXECUTABLE xsdcxx)

if (XSD_EXECUTABLE)
    file(GLOB_RECURSE XSD
            "${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader/*.xsd"
    )
    add_custom_target(xsd
            COMMAND ${XSD_EXECUTABLE} cxx-tree --std c++11 --hxx-suffix .h --cxx-suffix .cpp --generate-doxygen --generate-serialization ${XSD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader
            COMMENT "Generating code from xml schemas"
            VERBATIM
    )

endif()