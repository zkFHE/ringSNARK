#include <gtest/gtest.h>

#include "polynomials.hpp"
#include "seal/seal.h"
#include "../seal/seal_ring.hpp"
#include "evaluation_domain.hpp"
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
    class InterpolationTest : public testing::Test {

    };

    using Types = ::testing::Types<PrimitiveWrapper<double>, ringsnark::seal::RingElem>;
    TYPED_TEST_SUITE(InterpolationTest, Types);

    TYPED_TEST(InterpolationTest, TestInterpolation) {
        size_t n = 8;
        vector<TypeParam> coeffs(n);
        vector<TypeParam> x(n);
        vector<TypeParam> y(n);
        for (size_t i = 0; i < n; i++) {
            coeffs[i] = (i == 0) ? TypeParam::zero() : coeffs[i - 1] + TypeParam::one();
            x[i] = TypeParam(i);
        }
        for (size_t i = 0; i < n; i++) {
            y[i] = eval(coeffs, x[i]);
        }

        vector<TypeParam> coeffs_interpolated = interpolate<TypeParam>(x, y);
        vector<TypeParam> y_interpolated(n);
        for (size_t i = 0; i < n; i++) {
            y_interpolated[i] = eval(coeffs_interpolated, x[i]);
        }

        for (int i = 0; i < n; i++) {
            EXPECT_EQ(y[i], y_interpolated[i]);
        }

        for (int i = 0; i < n; i++) {
            EXPECT_EQ(coeffs[i], coeffs_interpolated[i]);
        }
    }

    TYPED_TEST(InterpolationTest, TestLagrangePolynomials) {
        size_t n = 8;
        auto domain = ringsnark::get_evaluation_domain<TypeParam>(n);
        vector<TypeParam> coeffs(n);
        vector<TypeParam> x(n);
        vector<TypeParam> y(n);
        for (size_t i = 0; i < n; i++) {
            coeffs[i] = (i == 0) ? TypeParam::zero() : coeffs[i - 1] + TypeParam::one();
            x[i] = domain->get_domain_element(i);
        }
        for (size_t i = 0; i < n; i++) {
            y[i] = eval(coeffs, x[i]);
        }

        for (size_t i = 0; i < 20; i++) {
            TypeParam s(domain->m + i);
            vector<TypeParam> lagrange = domain->evaluate_all_lagrange_polynomials(s);

            TypeParam res_interpolated = TypeParam::zero();
            for (size_t j = 0; j < n; j++) {
                res_interpolated += y[j] * lagrange[j];
            }
            TypeParam res = eval(coeffs, s);

            EXPECT_EQ(res_interpolated, res);
        }
    }

}

int main(int argc, char **argv) {
    ::seal::SEALContext context = get_context();
    ringsnark::seal::RingElem::set_context(context);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}