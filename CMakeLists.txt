cmake_minimum_required(VERSION 3.10)

# Define compilaton mode
# ON: release mode, OFF: debug model
option(MODE_release "release" ON) 

IF(MODE_release)
        set(CMAKE_CXX_FLAGS "-O3")
        set(FLUXOS_EXEC_name "fluxos")
ELSE()
        set(CMAKE_CXX_FLAGS "-g3 -Wall -pedantic")
        set(FLUXOS_EXEC_name "fluxos_debug")
ENDIF()

IF (UNIX)
	IF (UNIX AND NOT APPLE)
		MESSAGE(STATUS "current platform : Linux")
	ELSEIF(APPLE)
		MESSAGE(STATUS "current platform : MacOS")
	ENDIF()

        project(fluxos)

        file(GLOB SOURCES "src/*.cpp")

        set(CMAKE_CXX_STANDARD 17)
        
        # Packages
        find_package(Armadillo REQUIRED)
        find_package(OpenMP)


        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

        include_directories(${AMRADILLO_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})

        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

        # Create Results folder if non-existant
        file(MAKE_DIRECTORY Results)
        
        add_executable(${FLUXOS_EXEC_name}
                ${SOURCES})

        target_link_libraries(
                ${FLUXOS_EXEC_name}
                ${ARMADILLO_LIBRARIES})
        
	
ELSE()	
	MESSAGE(STATUS "curent platform : windows")
	project(fluxos)
	file(GLOB SOURCES "src/*.cpp")
	# Create Results folder if non-existant
	file(MAKE_DIRECTORY Results)

	set(CMAKE_CXX_STANDARD 17)

	# Packages
	find_package(OpenMP)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Results)
	add_executable(${FLUXOS_EXEC_name}, ${SOURCES})
	target_include_directories(${FLUXOS_EXEC_name} PUBLIC "I:/FLUXOS_cpp-master/Armadillo/include")

ENDIF()
