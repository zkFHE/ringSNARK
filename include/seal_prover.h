#ifndef RINOCCHIO_SEAL_PROVER_H
#define RINOCCHIO_SEAL_PROVER_H

#include "seal/seal.h"
#include "ring.h"
#include "ring_snark.h"
#include "poly_arith.h"

using namespace rinocchio;

class SealProver : public Prover<SealPoly, SealPoly, uint64_t> {
public:
    SealProver() = default;

protected:
    SealPoly multiply_inplace(SealPoly &e, const uint64_t &a) override;

    SealPoly multiply_inplace(SealPoly &e, const SealPoly &r) override;

    SealPoly add_inplace(SealPoly &e1, const SealPoly &e2) override;
};

#endif //RINOCCHIO_SEAL_PROVER_H
