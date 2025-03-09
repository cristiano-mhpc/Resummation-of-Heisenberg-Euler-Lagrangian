cmake_minimum_required(VERSION 3.14)
project(AutoLU LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find MPI
find_package(MPI REQUIRED)

# Find Boost with required components
find_package(Boost REQUIRED COMPONENTS mpi serialization thread system chrono)

# Find MPFR
find_path(MPFR_INCLUDE_DIR mpfr.h)
find_library(MPFR_LIBRARIES NAMES mpfr)

# Find GMP
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARIES NAMES gmp)

# Check if found
if (NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARIES)
    message(FATAL_ERROR "MPFR not found. Please install it and set MPFR_INCLUDE_DIR/MPFR_LIBRARIES.")
endif()
if (NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
    message(FATAL_ERROR "GMP not found. Please install it and set GMP_INCLUDE_DIR/GMP_LIBRARIES.")
endif()

# Include directories
include_directories(${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})

# Define the executable
add_executable(auto LU.cpp)

# Link libraries
target_link_libraries(auto PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
    pthread
)

# Optimization flags
target_compile_options(auto PRIVATE -O3)

