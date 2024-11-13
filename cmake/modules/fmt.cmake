set(FMT_VERSION 11.0.2)

if(NOT TARGET fmt)
    message("Fetching fmt : version 11.0.2")
    include(FetchContent)
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG        11.0.2
    )
    FetchContent_MakeAvailable(fmt)
endif()