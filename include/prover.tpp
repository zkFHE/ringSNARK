using namespace std;
using namespace rinocchio;

template<typename E, typename R, typename A>
Proof<E> Prover<E, R, A>::prove(SnarkParameters<R, A> params, CRS<E> crs) {
    size_t k = params.indices_mid[0];
    size_t i = 0;

    // Init A, A'
    R coeff = R(params.values[k]);
    coeff = multiply_inplace(coeff, params.v[k][i]);

    E a(crs.s_pows[i]);
    a = multiply_inplace(a, coeff);

    E a_alpha(crs.alpha_s_pows[i]);
    a_alpha = multiply_inplace(a_alpha, coeff);

    // B, B'
    coeff = R(params.values[k]);
    coeff = multiply_inplace(coeff, params.w[k][i]);

    E b(crs.s_pows[i]);
    b = multiply_inplace(b, coeff);

    E b_alpha(crs.alpha_s_pows[i]);
    b_alpha = multiply_inplace(b_alpha, coeff);

    // C, C'
    coeff = R(params.values[k]);
    coeff = multiply_inplace(coeff, params.y[k][i]);

    E c(crs.s_pows[i]);
    c = multiply_inplace(c, coeff);

    E c_alpha(crs.alpha_s_pows[i]);
    c_alpha = multiply_inplace(c_alpha, coeff);

    for (size_t i_k = 0; i_k < params.indices_mid.size(); i_k++) {
        size_t deg = params.v[k].size();
        k = params.indices_mid[i_k];
        for (i = (i_k == 0) ? 1 : 0; i < deg; i++) {
            // A, A'
            coeff = R(params.values[k]);
            coeff = multiply_inplace(coeff, params.v[k][i]);

            E tmp(crs.s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(a, tmp);

            tmp = E(crs.alpha_s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(a_alpha, tmp);

            // B, B'
            coeff = R(params.values[k]);
            coeff = multiply_inplace(coeff, params.w[k][i]);

            tmp = E(crs.s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(b, tmp);

            tmp = E(crs.alpha_s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(b_alpha, tmp);

            // C, C'
            coeff = R(params.values[k]);
            coeff = multiply_inplace(coeff, params.y[k][i]);

            tmp = E(crs.s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(c, tmp);

            tmp = E(crs.alpha_s_pows[i]);
            tmp = multiply_inplace(tmp, coeff);
            add_inplace(c_alpha, tmp);
        }
    }

    // Compute D, E
    E d = inner_prod_enc_ring(params.h.coeffs, crs.s_pows);
    E d_alpha = inner_prod_enc_ring(params.h.coeffs, crs.alpha_s_pows);

    // Compute F
    E f = inner_prod_enc_ring(params.values_mid, crs.beta_prod);

    return Proof<R>{a, a_alpha, b, b_alpha, c, c_alpha, d, d_alpha, f};
}

template<typename E, typename R, typename A>
E Prover<E, R, A>::inner_prod_enc_ring(const vector<R> &poly, const vector<E> &values) {
    assert(poly.size() <= values.size()); // poly might have lower degree
    E res(values[0]);
    res = multiply_inplace(res, poly[0]);
    for (size_t j = 1; j < poly.size(); j++) {
        // TODO: add check for poly == 0 \in R here?
        E tmp(values[j]);
        tmp = multiply_inplace(tmp, poly[j]);
        res = add_inplace(res, tmp);
    }
    return res;
}
