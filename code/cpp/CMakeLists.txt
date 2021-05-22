project(project_1_fys4411)

# Set CXX standard to C++11
set(CMAKE_CXX_STANDARD 11)

# Check if the user has specified a build type via the command line, e.g., (in the
# build-directory)
#
set(CMAKE_BUILD_TYPE Release)
#
# if not, set Debug as the default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Add all compiler warnings. Remember, a warning in C++ can very often lead to a bug
# later on so _fix your warnings!_
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# Add debugging flag for tools such as gdb and valgrind.
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add full compiler optimization when in the release build type. Other options include
# -O0, -O1, -O2.
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add directories containing header-files
include_directories(
    .
    Hamiltonians
    InitialStates
    Math
    WaveFunctions
)

# Add sources from all directories
file(
    GLOB SOURCES "*.cpp"
    GLOB SOURCES "Hamiltonians/*.cpp"
    GLOB SOURCES "InitialStates/*.cpp"
    GLOB SOURCES "Math/*.cpp"
    GLOB SOURCES "WaveFunctions/*.cpp"
)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})



add_library(mylib STATIC ${SOURCES})
find_package(OpenMP REQUIRED)



# Create executables

add_executable(vmc main/main.cpp)
add_executable(importance_sampling main/importance_sampling.cpp)
add_executable(coulomb_main main/coulomb_main.cpp)
add_executable(coulomb_importance_sampling main/coulomb_importance_sampling.cpp)

# Link libraries

target_link_libraries(vmc mylib)
target_link_libraries(vmc ${ARMADILLO_LIBRARIES})
target_link_libraries(vmc OpenMP::OpenMP_CXX)

target_link_libraries(importance_sampling mylib)
target_link_libraries(importance_sampling ${ARMADILLO_LIBRARIES})

target_link_libraries(coulomb_main mylib)
target_link_libraries(coulomb_main ${ARMADILLO_LIBRARIES})

target_link_libraries(coulomb_importance_sampling mylib)
target_link_libraries(coulomb_importance_sampling ${ARMADILLO_LIBRARIES})
