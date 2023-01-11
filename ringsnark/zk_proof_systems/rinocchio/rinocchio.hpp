#ifndef RINOCCHIO_HPP
#define RINOCCHIO_HPP

#include <ringsnark/zk_proof_systems/r1cs_ppzksnark.hpp>

using std::vector;

namespace ringsnark::rinocchio {
    template<typename RingT, typename EncT>
    class proving_key : public ringsnark::proving_key<RingT, EncT> {
        using PublicKey = typename EncT::PublicKey;

    public:
        const r1cs_constraint_system<RingT> constraint_system;
        const vector<EncT> s_pows, alpha_s_pows;
        const vector<EncT> beta_prods;
        const EncT beta_rv_ts, beta_rw_ts, beta_ry_ts;
        const EncT alpha_rv_ts, alpha_rw_ts, alpha_ry_ts; // Unused?
        const vector<EncT> rv_vs, rw_ws, ry_ys; // Unused?
        const PublicKey pk_enc;

        proving_key(const ringsnark::r1cs_constraint_system<RingT> &constraint_system,
                    const vector<EncT> &s_pows,
                    const vector<EncT> &alpha_s_pows,
                    const vector<EncT> &beta_prods,
                    const EncT beta_rv_ts, const EncT beta_rw_ts, const EncT beta_ry_ts,
                    const EncT alpha_rv_ts, const EncT alpha_rw_ts, const EncT alpha_ry_ts,
                    const vector<EncT> rv_vs, const vector<EncT> rw_ws, const vector<EncT> ry_ys,
                    const PublicKey &pk_enc) :
                constraint_system(constraint_system),
                s_pows(s_pows),
                alpha_s_pows(alpha_s_pows),
                beta_prods(beta_prods),
                beta_rv_ts(beta_rv_ts), beta_rw_ts(beta_rw_ts), beta_ry_ts(beta_ry_ts),
                alpha_rv_ts(alpha_rv_ts), alpha_rw_ts(alpha_rw_ts), alpha_ry_ts(alpha_ry_ts),
                rv_vs(rv_vs), rw_ws(rw_ws), ry_ys(ry_ys),
                pk_enc(pk_enc) {
            assert(s_pows.size() == constraint_system.num_constraints() + 1);
            assert(alpha_s_pows.size() == constraint_system.num_constraints() + 1);
            assert(beta_prods.size() == constraint_system.auxiliary_input_size);
        }

        [[nodiscard]] size_t size_in_bits() const override {
            return s_pows.size() * s_pows[0].size_in_bits()
                   + alpha_s_pows.size() * alpha_s_pows[0].size_in_bits()
                   + beta_prods.size() * beta_prods[0].size_in_bits()
                   + beta_rv_ts.size_in_bits() + beta_rw_ts.size_in_bits() + beta_ry_ts.size_in_bits()
                   + EncT::size_in_bits_pk(pk_enc);
        }
    };

    template<typename RingT, typename EncT>
    class verification_key : public ringsnark::verification_key<RingT, EncT> {
        using SecretKey = typename EncT::SecretKey;

    public:
        const ringsnark::rinocchio::proving_key<RingT, EncT> pk;
        const RingT s;
        const RingT alpha;
        const RingT beta;
        const RingT r_v, r_w, r_y;
        const SecretKey sk_enc{};

        verification_key() = default;

        verification_key(const proving_key<RingT, EncT> &pk, RingT s, RingT alpha, RingT beta,
                         RingT r_v, RingT r_w, RingT r_y,
                         SecretKey sk_enc) : pk(pk),
                                             s(s),
                                             alpha(alpha),
                                             beta(beta),
                                             r_v(r_v), r_w(r_w), r_y(r_y),
                                             sk_enc(sk_enc) {
            assert(alpha.is_invertible());
            assert(!beta.is_zero());
            assert(r_v.is_invertible());
            assert(r_w.is_invertible());
            assert(r_y == r_v * r_w);
        }

        [[nodiscard]] size_t size_in_bits() const override {
            return s.size_in_bits() + alpha.size_in_bits() + beta.size_in_bits() +
                   r_v.size_in_bits() + r_w.size_in_bits() + r_y.size_in_bits() +
                   EncT::size_in_bits_sk(sk_enc);
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
        const EncT A;
        const EncT A_prime;
        const EncT B;
        const EncT B_prime;
        const EncT C;
        const EncT C_prime;
        const EncT D;
        const EncT D_prime;
        const EncT F;

        proof(EncT a, EncT aPrime, EncT b, EncT bPrime, EncT c, EncT cPrime, EncT d, EncT dPrime, EncT f) : A(a),
                                                                                                            A_prime(aPrime),
                                                                                                            B(b),
                                                                                                            B_prime(bPrime),
                                                                                                            C(c),
                                                                                                            C_prime(cPrime),
                                                                                                            D(d),
                                                                                                            D_prime(dPrime),
                                                                                                            F(f) {}

        [[nodiscard]] size_t size_in_bits() const override {
            return A.size_in_bits() + A_prime.size_in_bits() + B.size_in_bits() + B_prime.size_in_bits() +
                   C.size_in_bits() + C_prime.size_in_bits() + D.size_in_bits() + D_prime.size_in_bits() +
                   F.size_in_bits();
        }

    };
}  // namespace ringsnark

#include "rinocchio.tcc"

#endif