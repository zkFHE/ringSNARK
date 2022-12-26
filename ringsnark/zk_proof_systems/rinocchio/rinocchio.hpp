#ifndef RING_SNARK_H
#define RING_SNARK_H

#include <utility>
#include <set>
#include "seal/seal.h"
#include <ringsnark/zk_proof_systems/r1cs_ppzksnark.hpp>
#include <ringsnark/relations/arithmetic_programs/qrp/qrp.hpp>

namespace ringsnark::rinocchio {
    template<typename RingT, typename EncT>
    class proving_key : public ringsnark::proving_key<RingT, EncT> {
        using PublicKey = typename EncT::PublicKey;

    public:
        const r1cs_constraint_system<RingT> constraint_system;
        const vector<EncT> s_pows, alpha_s_pows;    // size = |I_mid|
        const vector<EncT> beta_prods;               // size = |I_mid|
        const PublicKey pk_enc;

        proving_key(const ringsnark::r1cs_constraint_system<RingT> &constraint_system,
                    const vector<EncT> &sPows,
                    const vector<EncT> &alphaSPows,
                    const vector<EncT> &betaProds,
                    const PublicKey &pkEnc) :
                constraint_system(constraint_system),
                s_pows(sPows),
                alpha_s_pows(alphaSPows),
                beta_prods(betaProds),
                pk_enc(pkEnc) {
            assert(this->constraint_system.primary_input_size == constraint_system.primary_input_size);
            assert(this->constraint_system.auxiliary_input_size == constraint_system.auxiliary_input_size);
            assert(this->constraint_system.constraints == constraint_system.constraints);

//            assert(this->constraint_system == constraint_system);
        }

        [[nodiscard]] size_t size_in_bits() const override {
            return s_pows.size() * s_pows[0].size_in_bits()
                   + alpha_s_pows.size() * alpha_s_pows[0].size_in_bits()
                   + beta_prods.size() * beta_prods[0].size_in_bits();
            //TODO: + pk_enc.size_in_bits();
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

        verification_key(const proving_key<RingT, EncT> &pk, RingT s, RingT alpha, RingT beta, RingT rV, RingT rW,
                         RingT rY, SecretKey skEnc) : pk(pk),
                                                      s(s),
                                                      alpha(alpha),
                                                      beta(beta),
                                                      r_v(rV), r_w(rW), r_y(rY),
                                                      sk_enc(skEnc) {}

        [[nodiscard]] size_t size_in_bits() const override {
            return s.size_in_bits() + alpha.size_in_bits() + beta.size_in_bits() +
                   r_v.size_in_bits() + r_w.size_in_bits() + r_y.size_in_bits();
            // TODO: +sk_enc.size_in_bits();
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

//    template<typename RingT, typename EncT>
//    keypair<RingT, EncT> generator(const r1cs_constraint_system<RingT> &cs);
}  // namespace ringsnark

#include "rinocchio.tcc"

#endif