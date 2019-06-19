cmake_minimum_required(VERSION 3.10)
project(fluxos_cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-g3 -Wall -pedantic")
#set(CMAKE_CXX_FLAGS "-O2 -g ")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(Armadillo REQUIRED)
find_package(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${AMRADILLO_INCLUDE_DIR})

add_executable(fluxos_cpp
        main.cpp )

        #declavar.cpp)
target_link_libraries(
        fluxos_cpp
        ${ARMADILLO_LIBRARIES} )
