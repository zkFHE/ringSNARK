#ifndef RINGSNARK_SEAL_RING_HPP
#define RINGSNARK_SEAL_RING_HPP

#include <iostream>
#include <utility>
#include <vector>
#include <variant>
#include <memory>
#include "seal/seal.h"
#include "seal/util/rlwe.h"
#include "poly_arith.h"

using std::vector;

namespace ringsnark::seal {
    class RingElem {
    protected:
        using Scalar = uint64_t;
        using Poly = polytools::SealPoly;

        // TODO: try and turn this into a template parameter (which would probably require building a constexpr SEALContext)
        inline static ::seal::SEALContext *context = nullptr;
        std::variant<Poly, Scalar> value = (uint64_t) 0;
        inline static std::shared_ptr<::seal::UniformRandomGenerator> prng = nullptr;

        [[nodiscard]] Scalar &get_scalar();

    public:

        /*
         * Constructors
         */
        RingElem();

        RingElem(const RingElem &other) = default;

        RingElem(RingElem &&other) = default;

        RingElem &operator=(const RingElem &other) = default;

        virtual ~RingElem() = default;

        explicit RingElem(uint64_t value);

        explicit RingElem(const polytools::SealPoly &poly);

        /*
         * Static
         */
        static void set_context(::seal::SEALContext &context_) {
            if (context == nullptr) {
                context = new ::seal::SEALContext(context_);
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
            return RingElem(1);
        }

        static RingElem zero() {
            return RingElem(0);
        }

        static RingElem random_exceptional_element() {
            // TODO: throw error if number of exceptional elements is less than required
            auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();
            uint64_t q1 = parms.coeff_modulus()[0].value();
            // TODO: Switch to C++20 and use std::bit_width()
            uint64_t bit_width = 1ULL + (uint64_t) std::floor(std::log2l(q1));
            uint64_t mask = (1 << (bit_width + 1)) - 1;

            // Rejection sampling with masking
            uint64_t rand = ::seal::random_uint64() & mask;
            while (rand >= q1) { rand = ::seal::random_uint64() & mask; }

            return RingElem(rand);
        }

        inline static RingElem random_element() {
            if (prng == nullptr) {
                prng = ::seal::UniformRandomGeneratorFactory::DefaultFactory()->create();
            }

            auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();
            vector<uint64_t> coeffs(parms.poly_modulus_degree() * parms.coeff_modulus().size());
            ::seal::util::sample_poly_uniform(prng, parms, coeffs.data());
            return RingElem(polytools::SealPoly(get_context(), coeffs, get_context().first_parms_id()));
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

        [[nodiscard]] bool is_poly() const;

        [[nodiscard]] bool is_scalar() const;

        void negate_inplace();

        inline RingElem operator-() const {
            RingElem res(*this);
            res.negate_inplace();
            return res;
        }

        bool is_invertible() const noexcept;

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

        RingElem &to_poly_inplace();

        [[nodiscard]] RingElem &to_poly() const {
            auto *res = new RingElem(*this);
            res->to_poly_inplace();
            return *res;
        }

        [[nodiscard]] RingElem to_scalar() const;

        [[nodiscard]] size_t hash() const;

        class invalid_ring_elem_types : std::invalid_argument {
        public:
            explicit invalid_ring_elem_types() : invalid_argument("invalid types") {}
        };

        [[nodiscard]] Scalar get_scalar() const;

        [[nodiscard]] Poly get_poly() const;

        [[nodiscard]] Poly &get_poly();
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

        EncodingElem() = delete;

    public:
        using PublicKey = nullptr_t; // No keying material needed to evaluate affine combinations of ciphertexts
        using SecretKey = vector<::seal::SecretKey>;

        /*
         * Constructor
         */
        EncodingElem(const EncodingElem &other) : ciphertexts(other.ciphertexts) {
            assert(!other.ciphertexts.empty());
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

        static void set_context() {
            const ::seal::SEALContext &ring_context = RingElem::get_context();
            auto ring_params = ring_context.first_context_data()->parms();
            vector<::seal::SEALContext> enc_contexts;
            enc_contexts.reserve(ring_params.coeff_modulus().size());
            for (size_t i = 0; i < ring_params.coeff_modulus().size(); i++) {
                ::seal::EncryptionParameters enc_params(::seal::scheme_type::bgv);
                size_t poly_modulus_degree = ring_params.poly_modulus_degree();
                enc_params.set_poly_modulus_degree(poly_modulus_degree);
                enc_params.set_coeff_modulus(::seal::CoeffModulus::BFVDefault(poly_modulus_degree));
//                enc_params.set_plain_modulus(::seal::PlainModulus::Batching(poly_modulus_degree, ring_params));
                enc_params.set_plain_modulus(ring_params.plain_modulus());
                ::seal::SEALContext context(enc_params);

                if (context.first_context_data()->qualifiers().parameter_error !=
                    ::seal::EncryptionParameterQualifiers::error_type::success) {
                    std::cerr << context.first_context_data()->qualifiers().parameter_error_name() << std::endl;
                    std::cerr << context.first_context_data()->qualifiers().parameter_error_message() << std::endl;
                    throw std::invalid_argument("");
                }
                assert(context.first_context_data()->qualifiers().using_batching == true);
//                assert(context.using_keyswitching() == false); // TODO: can we always force this to be false while having enough noise budget for (potentially) many additions?

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

        static RingElem decode(const SecretKey &sk, const EncodingElem &e);

        /*
         * Members
         */
        [[nodiscard]] size_t size_in_bits() const;


        EncodingElem &operator+=(const EncodingElem &other);

        EncodingElem &operator*=(const RingElem &other);

        explicit EncodingElem(std::vector<::seal::Ciphertext> ciphertexts) : ciphertexts(std::move(ciphertexts)) {}
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
}

namespace std {
    template<>
    struct hash<ringsnark::seal::RingElem> {
        size_t operator()(const ringsnark::seal::RingElem &r) const {
            return r.hash();
        }
    };
}

#endif