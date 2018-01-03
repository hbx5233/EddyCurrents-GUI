cmake_minimum_required(VERSION 2.8)

project(EddyCurrents)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
if(EXISTS $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else()
  list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/v5.34.32/etc/cmake)
endif()
find_package(ROOT REQUIRED COMPONENTS RIO Net)

include_directories(${PROJECT_SOURCE_DIR}/headers
        ${ROOT_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/headers/*.h)

add_executable(EddyCurrents EddyCurrents.cpp ${sources} ${headers})
target_link_libraries(EddyCurrents ${ROOT_LIBRARIES})
