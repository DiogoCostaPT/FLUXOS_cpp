cmake_minimum_required(VERSION 3.10)

project(fluxos)

file(GLOB SOURCES "src/*.cpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-g3 -Wall -pedantic")
#set(CMAKE_CXX_FLAGS "-O2 -g ")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

# Packages
find_package(Armadillo REQUIRED)
find_package(OpenMP)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${AMRADILLO_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Create Results folder if non-existant
file(MAKE_DIRECTORY Results)

add_executable(fluxos
        ${SOURCES})

target_link_libraries(
        fluxos
        ${ARMADILLO_LIBRARIES})
