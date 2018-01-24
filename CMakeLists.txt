cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0054)
	cmake_policy(SET CMP0054 NEW)
endif()

project(STITCH C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions(-Wall -pedantic -fopenmp -O3)

find_package(Boost  REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")

add_subdirectory(BigWig)

add_subdirectory(core)

add_subdirectory(SPAN)

enable_testing()

add_subdirectory(test)
