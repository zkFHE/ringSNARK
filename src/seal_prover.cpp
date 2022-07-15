#include "seal_prover.h"

SealPoly SealProver::multiply_inplace(SealPoly &e, const uint64_t &a) {
    e.multiply_scalar_inplace(a);
    return e;
}

SealPoly SealProver::multiply_inplace(SealPoly &e, const SealPoly &r) {
    e.multiply_inplace(r);
    return e;
}

SealPoly SealProver::add_inplace(SealPoly& e1, const SealPoly& e2) {
    e1.add_inplace(e2);
    return e1;
}

