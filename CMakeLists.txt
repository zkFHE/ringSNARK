cmake_minimum_required(VERSION 3.20)

project(rinocchio LANGUAGES CXX)

# Library
find_package(SEAL 4.0 REQUIRED)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/depends")

add_library(rinocchio)
target_sources(rinocchio PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src/ring.cpp)
target_include_directories(rinocchio PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include/)


target_link_libraries(rinocchio PUBLIC polytools)
if (TARGET SEAL::seal)
    target_link_libraries(rinocchio PUBLIC SEAL::seal)
elseif (TARGET SEAL::seal_shared)
    target_link_libraries(rinocchio PUBLIC SEAL::seal_shared)
else ()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif ()

# Tests

# Example
add_executable(rinocchio-example)
target_sources(rinocchio-example PRIVATE ${CMAKE_CURRENT_LIST_DIR}/example.cpp)
target_link_libraries(rinocchio-example PRIVATE rinocchio)