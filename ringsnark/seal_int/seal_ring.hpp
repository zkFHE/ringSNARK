#ifndef RINGSNARK_SEAL_RING_INT_HPP
#define RINGSNARK_SEAL_RING_INT_HPP

#include "poly_arith.h"
#include "ringsnark/util/evaluation_domain.hpp"
#include "seal/seal.h"
#include "seal/util/uintarith.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintarithmod.h"
#include "seal/util/rlwe.h"
#include <iostream>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

using std::vector;

namespace ringsnark::seal_int {
    class RingElem {
    protected:
        // TODO: try and turn this into a template parameter (which would probably require building a constexpr SEALContext)
        inline static ::seal::SEALContext *context = nullptr;
        inline static std::shared_ptr<::seal::UniformRandomGenerator> prng = nullptr;
        friend std::ostream &operator<<(std::ostream &out, const RingElem &elem);
    public:
      std::vector<uint64_t> value;
      inline static vector<seal::Modulus> modulus;


        /*
         * Constructors
         */
        RingElem() = default;

        RingElem(const RingElem &other) = default;

        RingElem(RingElem &&other) = default;

        RingElem &operator=(const RingElem &other) = default;

        virtual ~RingElem() = default;

        RingElem(uint64_t x) {
          value = vector<uint64_t>(context->first_context_data()->parms().coeff_modulus().size());
          for (int i = 0; i < value.size(); i++) {
            value[i] = seal::util::barrett_reduce_64(x, context->first_context_data()->parms().coeff_modulus()[i]);
          }
        }

        RingElem(const vector<uint64_t>& value) : value(value) {
          assert(value.size() == context->first_context_data()->parms().coeff_modulus().size());
        }

        /*
         * Static
         */
        static void set_context(::seal::SEALContext &context_) {
            if (context == nullptr) {
                context = new ::seal::SEALContext(context_);
                modulus = context->first_context_data()->parms().coeff_modulus();
            } else {
                throw std::invalid_argument("cannot re-set context once set");
            }
        }

        static ::seal::SEALContext &get_context() {
            if (context == nullptr) {
                throw std::invalid_argument("context not set");
            } else {
                return *context;
            }
        }

        static RingElem one() {
            return RingElem(vector<uint64_t>(modulus.size(), 1));
        }

        static RingElem zero() {
            return RingElem(vector<uint64_t>(modulus.size(), 0));
        }

        static RingElem random_exceptional_element(const shared_ptr<evaluation_domain<RingElem>> domain = nullptr) {
            // TODO: throw error if number of exceptional elements is less than required
            auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();
            uint64_t q1 = parms.coeff_modulus()[0].value();
            uint64_t bit_width = 1ULL + (uint64_t) std::floor(std::log2l(q1));
            uint64_t mask = (1 << (bit_width + 1)) - 1;

            // Rejection sampling with masking
            uint64_t rand = ::seal::random_uint64() & mask;
            while (rand >= q1 || (domain && rand <= domain->m)) {
                rand = ::seal::random_uint64() & mask;
            }

            return RingElem(rand);
        }

        inline static RingElem random_element() {
            if (prng == nullptr) {
                prng = ::seal::UniformRandomGeneratorFactory::DefaultFactory()->create();
            }

            auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();
            vector<uint64_t> coeffs(parms.poly_modulus_degree() * parms.coeff_modulus().size());
            ::seal::util::sample_poly_uniform(prng, parms, coeffs.data());
            // TODO: FIXME
            return RingElem(polytools::SealPoly(get_context(), coeffs, &get_context().first_parms_id()).get_coefficient_rns(0));
        }

        static RingElem random_invertible_element() {
            RingElem res;
            do {
                res = random_element();
            } while (!res.is_invertible());
            return res;
        }

        static RingElem random_nonzero_element() {
            RingElem res;
            do {
                res = random_element();
            } while (res.is_zero());
            return res;
        }

        /*
         * Functions
         */
        [[nodiscard]] size_t size_in_bits() const;

        [[nodiscard]] bool is_zero() const;

        [[nodiscard]] bool fast_is_zero() const;

        void negate_inplace();

        inline RingElem operator-() const {
            RingElem res(*this);
            res.negate_inplace();
            return res;
        }

        [[nodiscard]] bool is_invertible() const noexcept;

        void invert_inplace();

        [[nodiscard]] inline RingElem inverse() const {
            RingElem res(*this);
            res.invert_inplace();
            return res;
        }

        RingElem &operator+=(const RingElem &other);

        RingElem &operator-=(const RingElem &other);

        RingElem &operator*=(const RingElem &other);

        RingElem &operator/=(const RingElem &other) {
            *this *= other.inverse();
            return *this;
        }

        [[nodiscard]] size_t hash() const;

        class invalid_ring_elem_types : std::invalid_argument {
        public:
            explicit invalid_ring_elem_types() : invalid_argument("invalid types") {}
        };
    };


    inline RingElem operator+(const RingElem &lhs, const RingElem &rhs) {
        RingElem res(lhs);
        res += rhs;
        return res;
    }

    inline RingElem operator-(const RingElem &lhs, const RingElem &rhs) {
        RingElem res(lhs);
        res -= rhs;
        return res;
    }

    inline RingElem operator*(const RingElem &lhs, const RingElem &rhs) {
        RingElem res(lhs);
        res *= rhs;
        return res;
    }

    inline RingElem operator/(const RingElem &lhs, const RingElem &rhs) {
        RingElem res(lhs);
        res /= rhs;
        return res;
    }

    bool operator==(const RingElem &lhs, const RingElem &rhs);

    inline bool operator!=(const RingElem &lhs, const RingElem &rhs) {
        return !operator==(lhs, rhs);
    }

    std::ostream &operator<<(std::ostream &out, const RingElem &elem);

    class EncodingElem {
    protected:
        inline static std::vector<::seal::SEALContext> contexts;
        // Use pointers instead of values as an ugly hack.
        // std::vectors requires a copy-constructor to be available, whereas SEAL deletes them for BatchEncoder and Evaluator
        inline static std::vector<::seal::BatchEncoder *> encoders;
        inline static std::vector<::seal::Evaluator *> evaluators;

        std::vector<::seal::Ciphertext> ciphertexts;

//        EncodingElem() = delete;

    public:
        using PublicKey = nullptr_t; // No keying material needed to evaluate affine combinations of ciphertexts
        using SecretKey = vector<::seal::SecretKey>;
        class decoding_error : std::invalid_argument {
          public:
            explicit decoding_error() : invalid_argument("decoding error") {}
        };

        /*
         * Constructor
         */
        EncodingElem() = default;

        EncodingElem(const EncodingElem &other) : ciphertexts(other.ciphertexts) {
            assert(!other.ciphertexts.empty());
        }

        [[nodiscard]] bool is_empty() const {
            return this->ciphertexts.empty();
        }

        /*
         * Static
         */
        static std::tuple<PublicKey, SecretKey> keygen() {
            SecretKey sk;
            sk.reserve(get_contexts().size());
            for (const auto &context: get_contexts()) {
                ::seal::KeyGenerator keygen(context);
                sk.push_back(keygen.secret_key());
            }
            PublicKey pk = nullptr;

            return {nullptr, sk};
        }

        static void set_context(size_t N=0) {
            // TODO: find (joint) primes Q_1, ..., Q_L for encoding schemes s.t.
            // Q_1 > q_l, and Q, resp. L are just barely big enough to allow for a linear homomorphism
            const ::seal::SEALContext &ring_context = RingElem::get_context();
            auto ring_params = ring_context.first_context_data()->parms();
            vector<::seal::SEALContext> enc_contexts;
            enc_contexts.reserve(ring_params.coeff_modulus().size());

            // Automagically find suitable poly_modulus_degree and coeff_modulus
//            auto max_plain_modulus = ring_params.coeff_modulus()[ring_params.coeff_modulus().size() - 1].value();
            //TODO: binary search to find optimal_poly_modulus_degree; this would require knowing how many additions need to be performed
            auto poly_modulus_degree = (N > 0) ? N : ring_params.poly_modulus_degree();
            auto coeff_modulus_max_bit_count = ::seal::CoeffModulus::MaxBitCount(poly_modulus_degree);
            vector<int> coeff_modulus_bit_counts;
            // TODO: distribute smoothly instead?
            while (coeff_modulus_max_bit_count > 60) {
                coeff_modulus_bit_counts.push_back(60);
                coeff_modulus_max_bit_count -= 60;
            }
            coeff_modulus_bit_counts.push_back(coeff_modulus_max_bit_count);
            auto coeff_modulus = ::seal::CoeffModulus::Create(poly_modulus_degree, coeff_modulus_bit_counts);

            for (size_t i = 0; i < ring_params.coeff_modulus().size(); i++) {
                ::seal::EncryptionParameters enc_params(::seal::scheme_type::bgv);
                enc_params.set_poly_modulus_degree(poly_modulus_degree);
                enc_params.set_plain_modulus(ring_params.coeff_modulus()[i].value());
                enc_params.set_coeff_modulus(coeff_modulus);
                ::seal::SEALContext context(enc_params);

                if (context.first_context_data()->qualifiers().parameter_error !=
                    ::seal::EncryptionParameterQualifiers::error_type::success) {
                    std::cerr << context.first_context_data()->qualifiers().parameter_error_name() << std::endl;
                    std::cerr << context.first_context_data()->qualifiers().parameter_error_message() << std::endl;
                    throw std::invalid_argument("");
                }

                assert(context.first_context_data()->qualifiers().using_batching == true);

//                assert(context.using_keyswitching() == false); //TODO: can we always force this to be false while having enough noise budget for (potentially) many additions?

                enc_contexts.push_back(context);
            }
            set_contexts(enc_contexts);
        }

        static void set_contexts(const vector<::seal::SEALContext> &contexts_) {
            if (contexts.empty()) {
                contexts = vector<::seal::SEALContext>(contexts_);
                encoders = vector<::seal::BatchEncoder *>();
                evaluators = vector<::seal::Evaluator *>();
                for (const ::seal::SEALContext &c: contexts) {
                    encoders.push_back(new ::seal::BatchEncoder(c));
                    evaluators.push_back(new ::seal::Evaluator(c));
                }
            } else {
                throw std::invalid_argument("cannot re-set contexts once set");
            }
        }

        static std::vector<::seal::SEALContext> &get_contexts() {
            if (contexts.empty()) {
                throw std::invalid_argument("context not set");
            } else {
                return contexts;
            }
        }


        // Encode all elements in rs using the same BatchEncoder and Encryptor objects for efficiency
        static std::vector<EncodingElem> encode(const SecretKey &sk, const std::vector<RingElem> &rs);

        static EncodingElem inner_product(
            vector<EncodingElem>::const_iterator a_start,
            vector<EncodingElem>::const_iterator a_end,
            vector<RingElem>::const_iterator b_start,
            vector<RingElem>::const_iterator b_end
            );

        static RingElem decode(const SecretKey &sk, const EncodingElem &e);

        /*
         * Members
         */
        [[nodiscard]] size_t size_in_bits() const;

        [[nodiscard]] static size_t size_in_bits_pk(const PublicKey &pk) {
            return 0;
        }

        [[nodiscard]] static size_t size_in_bits_sk(const SecretKey &sk) {
            size_t size = 0;
            for (size_t i = 0; i < sk.size(); i++) {
                auto context = get_contexts()[i];
                for (const auto &q_i: context.first_context_data()->parms().coeff_modulus()) {
                    size += q_i.bit_count() * context.first_context_data()->parms().poly_modulus_degree();
                }
            }
            return size;
        }

        EncodingElem &operator+=(const EncodingElem &other);

        EncodingElem &operator*=(const RingElem &other);

  explicit EncodingElem(std::vector<::seal::Ciphertext> ciphertexts)
      : ciphertexts(std::move(ciphertexts)) {}

  EncodingElem& modswitch_inplace();

        friend bool operator==(const EncodingElem &lhs, const EncodingElem &rhs);
    };

    inline EncodingElem operator+(const EncodingElem &lhs, const EncodingElem &rhs) {
        EncodingElem res(lhs);
        res += rhs;
        return res;
    }

    inline EncodingElem operator*(const EncodingElem &lhs, const RingElem &rhs) {
        EncodingElem res(lhs);
        res *= rhs;
        return res;
    }

    inline EncodingElem operator*(const RingElem &lhs, const EncodingElem &rhs) {
        EncodingElem res(rhs);
        res *= lhs;
        return res;
    }

    inline bool operator==(const EncodingElem &lhs, const EncodingElem &rhs) {
        if (lhs.ciphertexts.size() != rhs.ciphertexts.size()) {
            return false;
        }
        for (size_t i = 0; i < lhs.ciphertexts.size(); i++) {
            assert(lhs.ciphertexts[i].size() == rhs.ciphertexts[i].size());
            for (size_t j = 0; j < lhs.ciphertexts[i].size(); j++) {
                ::polytools::SealPoly l(RingElem::get_context(), lhs.ciphertexts[i], j);
                ::polytools::SealPoly r(RingElem::get_context(), rhs.ciphertexts[i], j);
                if (!l.is_equal(r)) {
                    return false;
                }
            }
        }
        return true;
    }
} // namespace ringsnark::seal_int


namespace std {
    template<>
    struct hash<ringsnark::seal_int::RingElem> {
        size_t operator()(const ringsnark::seal_int::RingElem &r) const {
            return r.hash();
        }
    };
} // namespace std

#include "seal_ring.tcc"

#endif