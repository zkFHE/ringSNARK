#include <gtest/gtest.h>

#include "polynomials.hpp"
#include "seal/seal.h"
#include "../seal/seal_ring.hpp"

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
    class PrimitiveWrapper {
    public:
        T val;

        static PrimitiveWrapper zero() { return PrimitiveWrapper(0); }

        static PrimitiveWrapper one() { return PrimitiveWrapper(1); }

        explicit PrimitiveWrapper(T val = 0) : val(val) {}

        PrimitiveWrapper &operator-() const {
            return *(new PrimitiveWrapper(-this->val));
        }

        PrimitiveWrapper &operator+=(const PrimitiveWrapper &rhs) {
            val += rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator-=(const PrimitiveWrapper &rhs) {
            val -= rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator*=(const PrimitiveWrapper &rhs) {
            val *= rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator/=(const PrimitiveWrapper &rhs) {
            val /= rhs.val;
            return *this;
        }

        bool operator==(const PrimitiveWrapper<T> &rhs) const {
            return abs(val - rhs.val) < 1e-3; // Approximate equality for double
        }

        bool operator!=(const PrimitiveWrapper<T> &rhs) const {
            return !this->operator==(rhs);
        }

    };

    template<typename T>
    inline PrimitiveWrapper<T> operator+(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator-(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator*(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator/(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T>
    inline std::ostream &operator<<(std::ostream &out, PrimitiveWrapper<T> const &w) {
        return out << w.val;
    }




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

}

int main(int argc, char **argv) {
    ::seal::SEALContext context = get_context();
    ringsnark::seal::RingElem::set_context(context);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}