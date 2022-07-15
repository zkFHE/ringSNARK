cmake_minimum_required(VERSION 3.20)

project(rinocchio LANGUAGES CXX)

# Library
find_package(SEAL 4.0 REQUIRED)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/depends")

add_library(rinocchio INTERFACE)
target_sources(rinocchio INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/seal_prover.cpp)
target_sources(rinocchio INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/seal_verifier.cpp)
target_include_directories(rinocchio INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include/)


target_link_libraries(rinocchio INTERFACE polytools)
if (TARGET SEAL::seal)
    target_link_libraries(rinocchio INTERFACE SEAL::seal)
elseif (TARGET SEAL::seal_shared)
    target_link_libraries(rinocchio INTERFACE SEAL::seal_shared)
else ()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif ()

# Tests

# Example
add_executable(rinocchio-example)
target_sources(rinocchio-example PRIVATE ${CMAKE_CURRENT_LIST_DIR}/example.cpp)
target_link_libraries(rinocchio-example PRIVATE rinocchio)