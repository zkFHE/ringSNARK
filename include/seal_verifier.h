#ifndef RINOCCHIO_SEAL_VERIFIER_H
#define RINOCCHIO_SEAL_VERIFIER_H

#include "seal/seal.h"
#include "ring.h"
#include "ring_snark.h"
#include "poly_arith.h"

using namespace rinocchio;
using E = SealPoly;
using R = SealPoly;
using A = uint64_t;

class SealVerifier : public Verifier<E, R, A> {
public:
    seal::EncryptionParameters sealParams;

    SealVerifier(SnarkParameters<R, A> snarkParams, seal::EncryptionParameters sealParams) :
            Verifier<E, R, A>(std::move(snarkParams)),
            sealParams(std::move(sealParams)) {}

    R decode(void *sk, const E &x) override;

    E encode(void *pk, const R &x) override;

    E encode(void *pk, const A &x) override;

    R multiply_inplace(R &r, const A &a) override;

    R multiply_inplace(R &r1, const R &r2) override;

    A multiply_inplace(A &a1, const A &a2) override;

    A add_inplace(A &a1, const A &a2) override;

    R add_inplace(R &r1, const R &r2) override;

    R subtract_inplace(R &r1, const R &r2) override;

    bool is_unit(const R &x) override;

    bool is_zero(const A &a) override;

    bool is_zero(const R &r) override;

    bool is_equal(const A &a1, const A &a2) override;

    bool is_equal(const R &r1, const R &r2) override;
};

#endif //RINOCCHIO_SEAL_VERIFIER_H
