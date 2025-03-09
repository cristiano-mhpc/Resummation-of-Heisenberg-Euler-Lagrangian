cmake_minimum_required(VERSION 3.10)
project(ConstantSolver)

# Enable C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Find Eigen (Assuming it's installed in /usr/include/eigen3)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Link MPFR and GMP
find_library(MPFR_LIB mpfr REQUIRED)
find_library(GMP_LIB gmp REQUIRED)

# Add executable
add_executable(constant constant.cpp)

# Link libraries
target_link_libraries(constant ${MPFR_LIB} ${GMP_LIB})

