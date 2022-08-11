project(memikrr)
cmake_minimum_required(VERSION 3.2)
set(TARGET ${PROJECT_NAME} CACHE STRING "output executable name")

# OpenSim uses C++11 language features.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and hook up to OpenSim.
# ----------------------------
set(OPENSIM_INSTALL_DIR "/home/lql/opensim_install")
find_package(OpenSim REQUIRED PATHS "${OPENSIM_INSTALL_DIR}")

set(EIGEN_INSTALL_DIR "/home/lql/software/opensim_dependencies_install/eigen")
find_package(Eigen3 REQUIRED PATHS "${EIGEN_INSTALL_DIR}")

set(KDL_INSTALL_DIR "/home/lql/orocos_kdl_install")
find_package(orocos_kdl REQUIRED PATHS "${KDL_INSTALL_DIR}")

set(KDL_PARSER_INSTALL_DIR "/home/lql/kdl_parser_install")
find_package(kdl_parser REQUIRED PATHS "${KDL_PARSER_INSTALL_DIR}")

# Configure this project.
# -----------------------
include_directories(${OpenSim_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${orocos_kdl_INCLUDE_DIRS})
include_directories(${kdl_parser_INCLUDE_DIRS})

include_directories(./inc)
file(GLOB SOURCE_FILES ./src/*.cpp)

add_executable(${TARGET} ${SOURCE_FILES} ./main.cpp)

target_link_libraries(${TARGET} ${OpenSim_LIBRARIES} ${orocos_kdl_LIBRARIES} ${kdl_parser_LIBRARIES})
