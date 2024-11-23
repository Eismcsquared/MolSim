find_program(XSD_EXECUTABLE NAMES xsdcxx)

if (XSD_EXECUTABLE)
    file(GLOB_RECURSE XSD
            "${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader/model/*.xsd"
    )
    add_custom_target(xsd
            COMMAND ${XSD_EXECUTABLE} cxx-tree --hxx-suffix=.h --cxx-suffix=.cpp --std c++20 --generate-doxygen --generate-serialization ${XSD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader/model
            VERBATIM
    )
endif()