cmake_minimum_required(VERSION 3.10)

cmake_minimum_required(VERSION 3.10)

project(GravitySolver C)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Add compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3 -g")

# Add MPI include directories
include_directories(${MPI_C_INCLUDE_DIRS})

# Add source files
set(SOURCE_FILES
    main.c
    grid.c
    tree.c
    gravity.c
)

# Add executable
add_executable(gravity_solver ${SOURCE_FILES})

# Link math library and MPI libraries
target_link_libraries(gravity_solver m ${MPI_LIBRARIES})
