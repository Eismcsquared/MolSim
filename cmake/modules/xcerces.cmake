# Check if xerces-c is already included
if (NOT TARGET xerces-c)
    
    include(FetchContent)
    FetchContent_Declare(
        xerces-c
        GIT_REPOSITORY https://github.com/apache/xerces-c.git
        GIT_TAG        v3.3.0 # v3.3.0 is the latest version at 07.11.2024
    )
    FetchContent_MakeAvailable(xerces-c)

    set(buildSamples OFF CACHE BOOL "Disable Xerces sample builds")
    set(buildTests OFF CACHE BOOL "Disable Xerces test builds")
endif()
