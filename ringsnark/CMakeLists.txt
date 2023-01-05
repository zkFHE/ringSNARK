add_library(ringsnark INTERFACE reductions/r1cs_to_qrp/evaluation_domain.tcc)

target_include_directories(ringsnark INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include/)

install(
        DIRECTORY "" DESTINATION "include/ringsnark"
        FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.tcc"
        PATTERN "*.cpp"
        PATTERN "tests" EXCLUDE
        PATTERN "examples" EXCLUDE
)

target_sources(ringsnark INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}/seal/seal_ring.hpp
        ${CMAKE_CURRENT_LIST_DIR}/seal/seal_ring.cpp
        ${CMAKE_CURRENT_LIST_DIR}/util/polynomials.hpp
        ${CMAKE_CURRENT_LIST_DIR}/util/polynomials.tcc
        ${CMAKE_CURRENT_LIST_DIR}/relations/variable.hpp
        ${CMAKE_CURRENT_LIST_DIR}/relations/variable.tcc
        ${CMAKE_CURRENT_LIST_DIR}/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp
        ${CMAKE_CURRENT_LIST_DIR}/relations/constraint_satisfaction_problems/r1cs/r1cs.tcc
        ${CMAKE_CURRENT_LIST_DIR}/relations/arithmetic_programs/qrp/qrp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/relations/arithmetic_programs/qrp/qrp.tcc
        ${CMAKE_CURRENT_LIST_DIR}/reductions/r1cs_to_qrp/evaluation_domain.hpp
        ${CMAKE_CURRENT_LIST_DIR}/reductions/r1cs_to_qrp/evaluation_domain.tcc
        ${CMAKE_CURRENT_LIST_DIR}/reductions/r1cs_to_qrp/r1cs_to_qrp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/reductions/r1cs_to_qrp/r1cs_to_qrp.tcc
        ${CMAKE_CURRENT_LIST_DIR}/gadgetlib/protoboard.hpp
        ${CMAKE_CURRENT_LIST_DIR}/gadgetlib/protoboard.tcc
        ${CMAKE_CURRENT_LIST_DIR}/zk_proof_systems/r1cs_ppzksnark.hpp
        )


add_executable(
        interpolation_test

        util/interpolation_test.cpp
        util/polynomials.hpp
        util/polynomials.tcc
        util/test_utils.hpp)

target_link_libraries(
        interpolation_test

        ringsnark
        gtest
        gtest_main
)

add_executable(
        division_test

        util/division_test.cpp
        util/polynomials.hpp
        util/polynomials.tcc
        util/test_utils.hpp)

target_link_libraries(
        division_test

        ringsnark
        gtest
        gtest_main
)


