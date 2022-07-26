#include "seal_prover.h"

R SealProver::multiply_inplace(R &r, const A &a) {
    r.multiply_scalar_inplace(a);
    return r;
}

E SealProver::multiply_inplace(E &e, const R &r) {
    e.multiply_inplace(r);
    return e;
}

E SealProver::add_inplace(E& e1, const E& e2) {
    e1.add_inplace(e2);
    return e1;
}

