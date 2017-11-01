cmake_minimum_required(VERSION 2.8)
project(carpio_1d_poisson)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/)

set(SOURCE_FILES main.cpp)

set(PROJECT_SOURCE_DIR ./)

set(DIR_CARPIO ../../..)

# location of header files
include_directories(
    ${DIR_CARPIO}/lib
)

add_subdirectory(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/utility 
    ${CMAKE_CURRENT_BINARY_DIR}/lib/utility
)

add_executable(main ${SOURCE_FILES})

# link the libraries to the executable
target_link_libraries (main utility)
