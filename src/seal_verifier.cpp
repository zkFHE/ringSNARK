#include "seal_verifier.h"

R SealVerifier::decode(void *sk, const E &x) {
    return x;
}

E SealVerifier::encode(void *pk, const R &x) {
    return x;
}

E SealVerifier::encode(void *pk, const A &x) {
    auto *context = static_cast<SEALContext *>(pk);
    Plaintext
            plain(seal::util::uint_to_hex_string(&x, 1));
    E res = E(*context, plain, &(context->first_parms_id()));
    res.ntt_inplace(context->get_context_data(res.get_parms_id())->small_ntt_tables());
    return res;
}

R SealVerifier::multiply_inplace(R &r, const A &a) {
    r.multiply_scalar_inplace(a);
    return r;
}

R SealVerifier::multiply_inplace(R &r1, const R &r2) {
    r1.multiply_inplace(r2);
    return r1;
}


A SealVerifier::multiply_inplace(A &a1, const A &a2) {
    // TODO: use bigint
    A q1 = sealParams.coeff_modulus()[0].value();
    return (a1 * a2) % q1;
}


A SealVerifier::add_inplace(A &a1, const A &a2) {
    // TODO: use bigint
    A q1 = sealParams.coeff_modulus()[0].value();
    return (a1 + a2) % q1;
}

R SealVerifier::add_inplace(R &r1, const R &r2) {
    r1.add_inplace(r2);
    return r1;
}

R SealVerifier::subtract_inplace(R &r1, const R &r2) {
    r1.subtract_inplace(r2);
    return r1;
}

bool SealVerifier::is_unit(const R &x) {
    // TODO: needed when generating some CRS parameters in R*
    return true; // TODO: use Euclidean algo over R_q1, ..., R_qL to check if inverse exists
}

bool SealVerifier::is_zero(const A &a) {
    A q1 = sealParams.coeff_modulus()[0].value();
    return (a % q1) == 0;
}

bool SealVerifier::is_zero(const R &r) {
    return r.is_zero();
}

bool SealVerifier::is_equal(const A &a1, const A &a2) {
    A q1 = sealParams.coeff_modulus()[0].value();
    return (a1 % q1) == (a2 % q1);
}

bool SealVerifier::is_equal(const R &r1, const R &r2) {
    return r1.is_equal(r2);
}
