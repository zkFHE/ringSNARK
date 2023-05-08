cmake_minimum_required(VERSION 3.20)

project(ringsnark LANGUAGES CXX)


add_executable(bench)
target_sources(bench PRIVATE ${CMAKE_CURRENT_LIST_DIR}/microbench.cpp)
target_link_libraries(bench PRIVATE ringsnark)

# Includes
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/depends")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ringsnark")

target_include_directories(ringsnark INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include/)

include_directories(.)

# Libraries
## lib: SEAL
# find_package(SEAL 4.0 REQUIRED)
if (TARGET SEAL::seal)
    target_link_libraries(ringsnark INTERFACE SEAL::seal)
elseif (TARGET SEAL::seal_shared)
    target_link_libraries(ringsnark INTERFACE SEAL::seal_shared)
else ()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif ()

find_package(Boost 1.80 COMPONENTS math)

## lib: polytools
target_link_libraries(ringsnark INTERFACE polytools)

## lib: Boost
target_include_directories(ringsnark INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(ringsnark INTERFACE ${Boost_LIBRARIES})

## lib: Google benchmark
target_link_libraries(ringsnark INTERFACE benchmark::benchmark)

# Executables
## Examples
add_executable(ringsnark_example_SEAL examples/example_SEAL.cpp)
target_link_libraries(ringsnark_example_SEAL PRIVATE ringsnark)

add_executable(ringsnark_example_plaintext_check_SEAL examples/example_plaintext_check_SEAL.cpp)
target_link_libraries(ringsnark_example_plaintext_check_SEAL PRIVATE ringsnark)

## Benchmarks
add_executable(ringsnark_bench_plaintext_check_SEAL benchmarks/bench_plaintext_check_SEAL.cpp)
target_link_libraries(ringsnark_bench_plaintext_check_SEAL PRIVATE ringsnark)
