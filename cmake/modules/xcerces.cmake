# Check if xerces-c is already included
if (NOT TARGET xerces-c)
    
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/_deps)

    FetchContent_Declare(
        xerces-c
        
        GIT_REPOSITORY https://github.com/apache/xerces-c.git
        GIT_TAG        v3.3.0 # v3.3.0 is the latest version at 07.11.2024
    )
    FetchContent_MakeAvailable(xerces-c)
endif()