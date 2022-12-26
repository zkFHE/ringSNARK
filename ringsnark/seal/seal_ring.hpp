#ifndef RINGSNARK_SEAL_RING_HPP
#define RINGSNARK_SEAL_RING_HPP

#include <iostream>
#include <utility>
#include <vector>
#include "seal/seal.h"
#include "poly_arith.h"

namespace ringsnark::seal {
    enum RING_ELEM_TYPE : size_t {
        SCALAR = 0,
        POLY = 1
    };

    class RingElem {
    protected:
        // TODO: try and turn this into a template parameter (which would probably require building a constexpr SEALContext)
        inline static ::seal::SEALContext *context = nullptr;
    public:
        polytools::SealPoly *poly = nullptr;
        uint64_t value = 0;
        RING_ELEM_TYPE type;

        /*
         * Constructors
         */
        RingElem();

        RingElem(const RingElem &other);

        RingElem(RingElem &&other) = default;

        RingElem &operator=(const RingElem &other) = default;

        RingElem &operator=(RingElem &&other) = default;

        virtual ~RingElem() = default;

        explicit RingElem(uint64_t value);

        RingElem(::seal::SEALContext &context_, ::seal::Plaintext &plaintext,
                 const ::seal::parms_id_type *parms_id_ptr);

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

        static RingElem random_invertible_element() {
            // TODO
            return RingElem::one();
        }

        static RingElem random_nonzero_element() {
            // TODO
            return RingElem::one();
        }

        /*
         * Functions
         */
        [[nodiscard]] size_t size_in_bits() const;

        bool is_zero() const;

        void negate_inplace();

        inline RingElem operator-() const {
            RingElem res(*this);
            res.negate_inplace();
            return res;
        }

        void invert_inplace();

        inline RingElem inverse() const {
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

        RingElem &to_poly() const {
            RingElem *res = new RingElem(*this);
            res->to_poly_inplace();
            return *res;
        }

        RingElem to_scalar() const;

        size_t hash() const;

        static void throw_invalid_types() {
            throw std::invalid_argument("invalid types");
        }
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
        inline static ::seal::SEALContext *context = nullptr;
        inline static ::seal::BatchEncoder *encoder = nullptr;
        inline static ::seal::Evaluator *evaluator = nullptr;

        std::vector<::seal::Ciphertext> ciphertexts;
    public:
        using PublicKey = nullptr_t; // No keying material needed to evaluate affine combinations of ciphertexts
        using SecretKey = ::seal::SecretKey;

        /*
         * Constructor
         */
        explicit EncodingElem(std::vector<::seal::Ciphertext> ciphertexts) : ciphertexts(std::move(ciphertexts)) {}

        /*
         * Static
         */
        static std::tuple<PublicKey, SecretKey> keygen() {
            ::seal::KeyGenerator keygen(get_context());
            const SecretKey &secret_key = keygen.secret_key();
            PublicKey pk = nullptr;
            return {nullptr, secret_key};
        }

        static void set_context(::seal::SEALContext &context_) {
            if (context == nullptr) {
                context = new ::seal::SEALContext(context_);
                encoder = new ::seal::BatchEncoder(context_);
                evaluator = new ::seal::Evaluator(context_);
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


        // Encode all elements in rs using the same BatchEncoder and Encryptor objects for efficiency
        static std::vector<EncodingElem> encode(const SecretKey &sk, const std::vector<RingElem> &rs);

        static RingElem decode(const SecretKey &sk, const EncodingElem& e);


        /*
         * Members
         */
        [[nodiscard]] size_t size_in_bits() const;


        EncodingElem &operator+=(const EncodingElem &other);

        EncodingElem &operator*=(const RingElem &other);

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