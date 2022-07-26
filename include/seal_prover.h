#ifndef RINOCCHIO_SEAL_PROVER_H
#define RINOCCHIO_SEAL_PROVER_H

#include "seal/seal.h"
#include "ring.h"
#include "ring_snark.h"
#include "poly_arith.h"

using namespace rinocchio;
using E = SealPoly;
using R = SealPoly;
using A = uint64_t;

class SealProver : public Prover<E, R, A> {
public:
    SealProver() = default;

protected:
    R multiply_inplace(R &r, const A &a) override;

    E multiply_inplace(E &e, const R &r) override;

    E add_inplace(E &e1, const E &e2) override;
};

#endif //RINOCCHIO_SEAL_PROVER_H
