# gtest.cmake

#check if gtest is already included
if (NOT TARGET gtest)
    # GoogleTest를 외부 프로젝트로 추가
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.15.2 # specify the version you need ,v1.15.2 is the latest version at 07.11.2021
    )
    FetchContent_MakeAvailable(googletest)

endif()
