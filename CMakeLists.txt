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


## lib: polytools
target_link_libraries(ringsnark INTERFACE polytools)

## lib: Boost
find_package(Boost 1.80 COMPONENTS math)
target_include_directories(ringsnark INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(ringsnark INTERFACE ${Boost_LIBRARIES})

## lib: Google benchmark
target_link_libraries(ringsnark INTERFACE benchmark::benchmark)

## lib: OpenFHE
set(CMAKE_CXX_STANDARD 17)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(OpenFHE)

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )

# Remove -Werror flag set by OpenFHE, since it prevents the compilation of other libs (e.g., Boost)
string(REPLACE " -Werror" " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()



# Executables
## Examples
add_executable(ringsnark_example_SEAL examples/example_SEAL.cpp)
target_link_libraries(ringsnark_example_SEAL PRIVATE ringsnark)

add_executable(ringsnark_example_OpenFHE examples/example_OpenFHE.cpp)
target_link_libraries(ringsnark_example_OpenFHE PRIVATE ringsnark)

add_executable(ringsnark_example_plaintext_check_SEAL examples/example_plaintext_check_SEAL.cpp)
target_link_libraries(ringsnark_example_plaintext_check_SEAL PRIVATE ringsnark)

## Benchmarks
add_executable(ringsnark_bench_plaintext_check_SEAL benchmarks/bench_plaintext_check_SEAL.cpp)
target_link_libraries(ringsnark_bench_plaintext_check_SEAL PRIVATE ringsnark)

add_executable(ringsnark_bench_mul_SEAL benchmarks/bench_mul_SEAL.cpp)
target_link_libraries(ringsnark_bench_mul_SEAL PRIVATE ringsnark)