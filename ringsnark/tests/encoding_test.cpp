#include <gtest/gtest.h>

#include "../seal/seal_ring.hpp"
#include "../util/test_utils.hpp"

size_t NUM_RING_ELEMS = 2;

::seal::SEALContext get_context() {
    ::seal::EncryptionParameters params(::seal::scheme_type::bgv);
    auto poly_modulus_degree = (size_t) pow(2, 11);
    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(::seal::CoeffModulus::BFVDefault(poly_modulus_degree));
    params.set_plain_modulus(::seal::PlainModulus::Batching(poly_modulus_degree, 20));
    ::seal::SEALContext context(params);
    return context;
}

namespace {
    template<typename T>
    class EncodingTest : public testing::Test {
    };

    using Types = ::testing::Types<std::pair<ringsnark::seal::RingElem, ringsnark::seal::EncodingElem>>;
    TYPED_TEST_SUITE(EncodingTest, Types);

    TYPED_TEST(EncodingTest, TestEncodingDecoding) {
        using RingT = typename TypeParam::first_type;
        using EncT = typename TypeParam::second_type;

        vector<RingT> rs;
        rs.reserve(NUM_RING_ELEMS);
        size_t i = 0;
        for (; i < NUM_RING_ELEMS / 2; i++) {
//            rs.push_back(RingT::random_exceptional_element());
            rs.push_back(RingT(5260053));
        }
        for (; i < NUM_RING_ELEMS; i++) {
            rs.push_back(RingT::random_element());
        }

        auto [pk, sk] = EncT::keygen();
        vector<EncT> es = EncT::encode(sk, rs);
        for (i = 0; i < NUM_RING_ELEMS; i++) {
            RingT r = EncT::decode(sk, es[i]);
            EXPECT_EQ(r, rs[i]);
        }
    }
}

int main(int argc, char **argv) {
    ::seal::SEALContext context = get_context();
    ringsnark::seal::RingElem::set_context(context);
    ringsnark::seal::EncodingElem::set_context();

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}