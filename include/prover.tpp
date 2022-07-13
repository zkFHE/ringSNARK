using namespace std;
using namespace rinocchio;

template<typename E, typename R, typename A>
Proof<E> Prover<E, R, A>::prove(SnarkParameters<R, A> params, CRS<E> crs) {
    // Evaluate V_k, V_k', ..., Y_k, Y_k'
    vector<E> v_s, w_s, y_s, v_alpha_s, w_alpha_s, y_alpha_s;
    v_s.reserve(params.num_mid);
    v_alpha_s.reserve(params.num_mid);
    w_s.reserve(params.num_mid);
    w_alpha_s.reserve(params.num_mid);
    y_s.reserve(params.num_mid);
    y_alpha_s.reserve(params.num_mid);
    for (const size_t k: params.indices_mid) {
        v_s.push_back(inner_prod_enc_exc(params.v[k], crs.s_pows));
        w_s.push_back(inner_prod_enc_exc(params.w[k], crs.s_pows));
        y_s.push_back(inner_prod_enc_exc(params.y[k], crs.s_pows));

        v_alpha_s.push_back(inner_prod_enc_exc(params.v[k], crs.alpha_s_pows));
        w_alpha_s.push_back(inner_prod_enc_exc(params.w[k], crs.alpha_s_pows));
        y_alpha_s.push_back(inner_prod_enc_exc(params.y[k], crs.alpha_s_pows));
    }

    // Evaluate A, A', ..., C, C'
    E a = inner_prod_enc_ring(params.values_mid, v_s);
    E a_alpha = inner_prod_enc_ring(params.values_mid, v_alpha_s);
    E b = inner_prod_enc_ring(params.values_mid, w_s);
    E b_alpha = inner_prod_enc_ring(params.values_mid, w_alpha_s);
    E c = inner_prod_enc_ring(params.values_mid, y_s);
    E c_alpha = inner_prod_enc_ring(params.values_mid, y_alpha_s);

    // Compute D, E
    E d = inner_prod_enc_ring(params.h.coeffs, crs.s_pows);
    E d_alpha = inner_prod_enc_ring(params.h.coeffs, crs.alpha_s_pows);

    // Compute F
    // TODO: do we need F for security?
    E f = inner_prod_enc_ring(params.values_mid, crs.beta_prod);

    return Proof<R>{a, a_alpha, b, b_alpha, c, c_alpha, d, d_alpha, f};
}

template<typename E, typename R, typename A>
E Prover<E, R, A>::inner_prod_enc_ring(const vector<R> &poly, const vector<E> &values) {
    assert(poly.size() <= values.size()); // poly might have lower degree
    E res(values[0]);
    res.multiply_inplace(poly[0]);
    for (size_t j = 1; j < poly.size(); j++) {
        // TODO: add check for poly == 0 \in R here?
        E tmp(values[j]);
        tmp.multiply_inplace(poly[j]);
        res.add_inplace(tmp);
    }
    return res;
}

template<typename E, typename R, typename A>
E Prover<E, R, A>::inner_prod_enc_exc(const vector<A> &poly, const vector<E> &values) {
    assert(poly.size() == values.size());
    E res(values[0]);
    res.multiply_scalar_inplace(poly[0]);
    for (size_t j = 1; j < values.size(); j++) {
        if (poly[j] != 0) {
            E tmp(values[j]);
            tmp.multiply_scalar_inplace(poly[j]);
            res.add_inplace(tmp);
        }
    }
    return res;
}
