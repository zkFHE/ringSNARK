#include <gtest/gtest.h>

#include "polynomials.hpp"
#include "seal/seal.h"
#include "../seal/seal_ring.hpp"
#include "test_utils.hpp"

::seal::SEALContext get_context() {
    ::seal::EncryptionParameters params(::seal::scheme_type::bgv);
    auto poly_modulus_degree = (size_t) pow(2, 10);
    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(::seal::CoeffModulus::BFVDefault(poly_modulus_degree));
    params.set_plain_modulus(::seal::PlainModulus::Batching(poly_modulus_degree, 20));
    ::seal::SEALContext context(params);
    return context;
}

namespace {
    template<typename T>
    class DivisionTest : public testing::Test {
    };

    using Types = ::testing::Types<PrimitiveWrapper<double>, ringsnark::seal::RingElem>;
    TYPED_TEST_SUITE(DivisionTest, Types);

    TYPED_TEST(DivisionTest, TestDivision) {
        size_t n = 110;

        vector<TypeParam> x(n);
        vector<TypeParam> quotient(n);
        for (size_t i = 0; i < n; i++) {
            x[i] = TypeParam(2 * i + 1); //TypeParam::random_element();
            quotient[i] = TypeParam(i + 1); //TypeParam::random_exceptional_element();
        }

        vector<TypeParam> y = multiply(quotient, x);

        vector<TypeParam> q = divide(y, x);

        EXPECT_LE(q.size(), quotient.size());
        while (q.size() < quotient.size()) {
            q.push_back(TypeParam::zero());
        }
        for (int i = 0; i < n; i++) {
            EXPECT_EQ(q[i], quotient[i]);
        }
    }
}

int main(int argc, char **argv) {
    ::seal::SEALContext context = get_context();
    ringsnark::seal::RingElem::set_context(context);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}