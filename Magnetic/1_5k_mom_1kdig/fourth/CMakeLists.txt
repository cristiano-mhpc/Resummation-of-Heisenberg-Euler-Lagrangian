cmake_minimum_required(VERSION 3.10)
project(FirstProject LANGUAGES CXX)

# Find required packages
find_package(MPI REQUIRED)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Boost
find_package(Boost REQUIRED COMPONENTS mpi serialization thread system chrono)

# Include directories
include_directories(${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

# Define the executable
add_executable(fourth   fourth.cpp  )

# Link libraries
target_link_libraries(fourth   PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    mpfr
    gmp
    pthread
)

# Optimization flags
target_compile_options(fourth   PRIVATE -Wall -Wextra -O2)

