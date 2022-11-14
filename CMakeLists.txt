cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(PROJECT_NAME hybridPNM)
project(${PROJECT_NAME})
add_subdirectory(../raplea ${PROJECT_BINARY_DIR}/raplea)

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../raplea)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../netgrid)

set(SOURCE_CODE Props.cpp Variables.cpp equation/Equation.cpp
        loc/Loc.cpp loc/LocP.cpp loc/LocS.cpp
        con/Con.cpp con/ConP.cpp con/ConS.cpp)

add_library(${PROJECT_NAME} ${SOURCE_CODE})
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen raplea)

pybind11_add_module(${PROJECT_NAME}_bind wrapper.cpp)
target_link_libraries(${PROJECT_NAME}_bind PRIVATE ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}_bind PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}_bind PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/..)

add_executable(${PROJECT_NAME}_demo ${SOURCE_CODE} demo/demo.cpp)
target_link_libraries(${PROJECT_NAME}_demo PUBLIC ${PROJECT_NAME})

