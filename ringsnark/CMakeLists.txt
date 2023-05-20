add_library(ringsnark INTERFACE util/evaluation_domain.tcc)

target_include_directories(ringsnark INTERFACE include/)

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
        seal_int/seal_ring.hpp
        seal_int/seal_ring.tcc
        seal/seal_ring.hpp
        seal/seal_ring.tcc
        seal/seal_util.hpp
        util/polynomials.hpp
        util/polynomials.tcc
        relations/variable.hpp
        relations/variable.tcc
        relations/constraint_satisfaction_problems/r1cs/r1cs.hpp
        relations/constraint_satisfaction_problems/r1cs/r1cs.tcc
        relations/arithmetic_programs/qrp/qrp.hpp
        relations/arithmetic_programs/qrp/qrp.tcc
        util/evaluation_domain.hpp
        util/evaluation_domain.tcc
        reductions/r1cs_to_qrp/r1cs_to_qrp.hpp
        reductions/r1cs_to_qrp/r1cs_to_qrp.tcc
        gadgetlib/protoboard.hpp
        gadgetlib/protoboard.tcc
        zk_proof_systems/r1cs_ppzksnark.hpp
        )


add_executable(
        interpolation_test

        util/interpolation_test.cpp
        util/polynomials.hpp
        util/polynomials.tcc
        util/test_utils.hpp
)

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
        util/test_utils.hpp
)

target_link_libraries(
        division_test

        ringsnark
        gtest
        gtest_main
)


add_executable(
        encoding_test

        tests/encoding_test.cpp
)

target_link_libraries(
        encoding_test

        ringsnark
        gtest
        gtest_main
)