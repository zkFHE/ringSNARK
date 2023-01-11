#ifndef GROTH16_HPP
#define GROTH16_HPP

#include <ringsnark/zk_proof_systems/r1cs_ppzksnark.hpp>

using std::vector;

namespace ringsnark::groth16 {
    template<typename RingT, typename EncT>
    class proving_key : public ringsnark::proving_key<RingT, EncT> {
        using PublicKey = typename EncT::PublicKey;

    public:
        const r1cs_constraint_system<RingT> constraint_system;
        const EncT alpha, beta;
        const vector<EncT> s_pows;
        const vector<EncT> gamma_io;
        const vector<EncT> delta_mid;
        const vector<EncT> delta_ts;
        const PublicKey pk_enc;

        proving_key(const ringsnark::r1cs_constraint_system<RingT> &constraint_system,
                    const EncT alpha,
                    const EncT beta,
                    const vector<EncT> &s_pows,
                    const vector<EncT> &gamma_io,
                    const vector<EncT> &delta_mid,
                    const vector<EncT> &delta_ts,
                    const PublicKey &pk_enc) :
                constraint_system(constraint_system),
                alpha(alpha),
                beta(beta),
                s_pows(s_pows),
                gamma_io(gamma_io),
                delta_mid(delta_mid),
                delta_ts(delta_ts),
                pk_enc(pk_enc) {
            assert(s_pows.size() == constraint_system.num_constraints() + 1);
            assert(gamma_io.size() == constraint_system.primary_input_size+1);
            assert(delta_mid.size() == constraint_system.auxiliary_input_size);
            assert(delta_ts.size() == constraint_system.num_constraints() + 1);
        }

        [[nodiscard]] size_t size_in_bits() const override {
            return s_pows.size() * s_pows[0].size_in_bits() +
                   gamma_io.size() * gamma_io[0].size_in_bits() +
                   delta_mid.size() * delta_mid[0].size_in_bits() +
                   delta_ts.size() * delta_ts[0].size_in_bits();
        }
    };

    template<typename RingT, typename EncT>
    class verification_key : public ringsnark::verification_key<RingT, EncT> {
        using SecretKey = typename EncT::SecretKey;

    public:
        const ringsnark::groth16::proving_key<RingT, EncT> pk;
        const RingT s;
        const RingT alpha;
        const RingT beta;
        const RingT gamma;
        const RingT delta;
        const SecretKey sk_enc{};

        verification_key() = default;

        verification_key(const proving_key<RingT, EncT> &pk, RingT s, RingT alpha, RingT beta, RingT gamma, RingT delta,
                         SecretKey sk_enc) : pk(pk),
                                             s(s),
                                             alpha(alpha),
                                             beta(beta),
                                             gamma(gamma),
                                             delta(delta),
                                             sk_enc(sk_enc) {
//            assert(s.is_exceptional());
            assert(alpha.is_invertible());
            assert(beta.is_invertible());
            assert(gamma.is_invertible());
            assert(delta.is_invertible());
        }

        [[nodiscard]] size_t size_in_bits() const override {
            return s.size_in_bits() + alpha.size_in_bits() + beta.size_in_bits() +
                   gamma.size_in_bits() + delta.size_in_bits() + EncT::size_in_bits_sk(sk_enc);
        }
    };

    template<typename RingT, typename EncT>
    class keypair {
    public:
        const proving_key<RingT, EncT> &pk;
        const verification_key<RingT, EncT> &vk;

        keypair(const keypair<RingT, EncT> &other) = default;

        keypair(keypair<RingT, EncT> &&other) noexcept = default;

        keypair(proving_key<RingT, EncT> &pk, verification_key<RingT, EncT> &vk) : pk(pk), vk(vk) {}

        keypair(proving_key<RingT, EncT> &&pk, verification_key<RingT, EncT> &&vk) : pk(std::move(pk)),
                                                                                     vk(std::move(vk)) {}
    };

    template<typename RingT, typename EncT>
    class proof : public ringsnark::proof<RingT, EncT> {
    public:
        const EncT A, B, C;

        proof(EncT a, EncT b, EncT c) : A(a), B(b), C(c) {}

        [[nodiscard]] size_t size_in_bits() const override {
            return A.size_in_bits() + B.size_in_bits() + C.size_in_bits();
        }

    };
}  // namespace ringsnark

#include "groth16.tcc"

#endif