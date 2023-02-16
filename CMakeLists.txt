cmake_minimum_required(VERSION 3.1)
project(multi_mc)

set(CMAKE_CXX_FLAGS "-Ofast -mtune=native -fdiagnostics-color=always")

# list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
# include(${ROOT_USE_FILE})

# include_directories(${ROOT_INCLUDE_DIRS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)

#add_library(mc_libs SHARED ${sources} ${headers})
#target_link_libraries(mc_libs ${ROOT_LIBRARIES})
add_executable(multi_mc main.cpp)
target_link_libraries(multi_mc ${ROOT_LIBRARIES})

#target_link_libraries(multi_mc mc_libs)