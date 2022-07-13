using namespace std;
using namespace rinocchio;

template<typename E, typename R, typename A>
bool Verifier<E, R, A>::is_unit(R x) {
    // TODO: needed when generating some CRS parameters in R*
    return true; // TODO: use Euclidean algo over R_q1, ..., R_qL to check if inverse exists
}

template<typename E, typename R, typename A>
VerifKey<E, R, A>
Verifier<E, R, A>::generate_vk(SnarkParameters<R, A> snark_parameters, A s, R alpha, R beta, R r_v, R r_w, R r_y) {
    // ({E(s^i)}_{i=0}^{num_mid}}, {E(alpha * s^i)}_{i=0}^{num_mid}}, {beta_prod}_{i=0}^{num_mid}, pk)
    vector<E> s_pows, alpha_s_pows, beta_prods;
    A curr_s = 1;
    R curr_alpha_s(alpha);
    assert(alpha.is_ntt_form());
    assert(beta.is_ntt_form());
    assert(r_v.is_ntt_form());
    assert(r_w.is_ntt_form());
    assert(r_y.is_ntt_form());

    assert(is_unit(alpha));
    assert(is_unit(r_v));
    assert(is_unit(r_w));
    assert(is_unit(r_y));

    assert(!beta.is_zero());

    uint64_t q1 = encryption_parameters.coeff_modulus()[0].value();
    for (size_t i = 0; i < snark_parameters.v.size(); i++) {
        auto is_in_Zq1 = [q1](uint64_t x) { return x < q1; };
        assert(all_of(snark_parameters.v[i].begin(), snark_parameters.v[i].end(), is_in_Zq1));
        assert(all_of(snark_parameters.w[i].begin(), snark_parameters.w[i].end(), is_in_Zq1));
        assert(all_of(snark_parameters.y[i].begin(), snark_parameters.y[i].end(), is_in_Zq1));
    }

    s_pows.push_back(encode(&snark_parameters.context, curr_s));
    alpha_s_pows.push_back(encode(&snark_parameters.context, curr_alpha_s));
    for (size_t k: snark_parameters.indices_mid) {
        curr_s = (curr_s * s) % q1;
        curr_alpha_s.multiply_scalar_inplace(s);
        s_pows.push_back(encode(&snark_parameters.context, curr_s));
        alpha_s_pows.push_back(encode(&snark_parameters.context, curr_alpha_s));

        A v_k_s = eval(snark_parameters.v[k], s);
        A w_k_s = eval(snark_parameters.w[k], s);
        A y_k_s = eval(snark_parameters.y[k], s);

        R beta_prod = R(r_v);
        beta_prod.multiply_scalar_inplace(v_k_s);

        R tmp = R(r_w);
        tmp.multiply_scalar_inplace(w_k_s);
        beta_prod.add_inplace(tmp);

        tmp = R(r_y);
        tmp.multiply_scalar_inplace(y_k_s);
        beta_prod.add_inplace(tmp);

        beta_prod.multiply_inplace(beta);

        beta_prods.push_back(encode(&snark_parameters.context, beta_prod));
    }
    CRS<E> crs{s_pows, alpha_s_pows, beta_prods};
    return VerifKey<E, R, A>{nullptr, crs, s, alpha, beta, r_v, r_w, r_y};
}

template<typename E, typename R, typename A>
bool Verifier<E, R, A>::verify(VerifKey<E, R, A> vk, SnarkParameters<R, A> params, Proof<E> proof) {
    // A = E(V_mid)
    // A' = E(V'_mid)
    // B = E(W_mid)
    // B' = E(W'_mid)
    // C = E(Y_mid)
    // C' = E(Y'_mid)
    // D = E(H)
    // D' = E(H')
    // F = E(L)
    R V_mid = decode(vk.sk, proof.A);
    R V_mid_prime = decode(vk.sk, proof.A_prime);
    R W_mid = decode(vk.sk, proof.B);
    R W_mid_prime = decode(vk.sk, proof.B_prime);
    R Y_mid = decode(vk.sk, proof.C);
    R Y_mid_prime = decode(vk.sk, proof.C_prime);
    R H = decode(vk.sk, proof.D);
    R H_prime = decode(vk.sk, proof.D_prime);
    R L_beta = decode(vk.sk, proof.F);

    // v_io(x) = v_0(x) + sum_{k=1}^l a_k * v_k(x)
    // w_io(x) = w_0(x) + sum_{k=1}^l a_k * w_k(x)
    // y_io(x) = y_0(x) + sum_{k=1}^l a_k * y_k(x)


    R *v_io_s = nullptr;
    R *w_io_s = nullptr;
    R *y_io_s = nullptr;
    for (const size_t k: params.indices_io) {
        const auto c_k = params.values[k];
        A v_k_s = eval(params.v[k], vk.s);
        if (v_k_s != 0) {
            R tmp = R(c_k);
            tmp.multiply_scalar_inplace(v_k_s);
            if (v_io_s) {
                v_io_s->add_inplace(tmp);
            } else {
                v_io_s = new R(tmp);
            }
        }

        A w_k_s = eval(params.w[k], vk.s);
        if (w_k_s != 0) {
            R tmp = R(c_k);
            tmp.multiply_scalar_inplace(w_k_s);
            if (w_io_s) {
                w_io_s->add_inplace(tmp);
            } else {
                w_io_s = new R(tmp);
            }
        }

        A y_k_s = eval(params.y[k], vk.s);
        if (y_k_s != 0) {
            R tmp = R(c_k);
            tmp.multiply_scalar_inplace(y_k_s);
            if (y_io_s) {
                y_io_s->add_inplace(tmp);
            } else {
                y_io_s = new R(tmp);
            }
        }
    }

    // L_span = (vk.r_v * V_mid) + (vk.r_w * W_mid) + (vk.r_y * Y_mid)
    R L = R(V_mid);
    L.multiply_inplace(vk.r_v);

    R tmp = (W_mid);
    tmp.multiply_inplace(vk.r_w);
    L.add_inplace(tmp);

    tmp = R(Y_mid);
    tmp.multiply_inplace(vk.r_y);
    L.add_inplace(tmp);

    // P = (v_io(s) + V_mid) * (w_io(s) + W_mid) - (y_io(s) + Y_mid)
    // R P = (v_io_s + V_mid) * (w_io_s + W_mid) - (y_io_s + Y_mid);
    R P = R(V_mid);
    if (v_io_s) {
        P.add_inplace(*v_io_s);
    }

    tmp = R(W_mid);
    if (w_io_s) {
        tmp.add_inplace(*w_io_s);
    }
    P.multiply_inplace(tmp);

    if (y_io_s) {
        P.subtract_inplace(*y_io_s);
    }
    P.subtract_inplace(Y_mid);

    // CHECK: V'_mid = alpha * V_mid
    tmp = R(V_mid);
    tmp.multiply_inplace(vk.alpha);
    if (!V_mid_prime.is_equal(tmp)) {
        return false;
    }
    // CHECK: W'_mid = alpha * W_mid
    tmp = R(W_mid);
    tmp.multiply_inplace(vk.alpha);
    if (!W_mid_prime.is_equal(tmp)) {
        return false;
    }
    // CHECK: Y'_mid = alpha * Y_mid
    tmp = R(Y_mid);
    tmp.multiply_inplace(vk.alpha);
    if (!Y_mid_prime.is_equal(tmp)) {
        return false;
    }
    // CHECK: H' = alpha * H
    tmp = R(H);
    tmp.multiply_inplace(vk.alpha);
    if (!H_prime.is_equal(tmp)) {
        return false;
    }
    // CHECK: L = beta * L_span
    tmp = R(L);
    tmp.multiply_inplace(vk.beta);
    if (!L_beta.is_equal(tmp)) {
        return false;
    }
    // P = H * t(s)
    R H_tmp = R(H);
    H_tmp.multiply_scalar_inplace(eval(params.t, vk.s));
    if (!P.is_equal(H_tmp)) {
        return false;
    }
    return true;
}

template<typename E, typename R, typename A>
A Verifier<E, R, A>::eval(const vector<A> &poly, const A x) {
    A q1 = encryption_parameters.coeff_modulus()[0].value();
    A res = poly[poly.size() - 1];
    for (int i = poly.size() - 2; i >= 0; i--) {
        res = (res * x) % q1;
        res = (res + poly[i]) % q1;
    }
    return res;
}


template<typename E, typename R, typename A>
E Verifier<E, R, A>::encode(void *pk, R x) {
//        SEALContext *context = static_cast<SEALContext *>(pk);
//        x.ntt_inplace(context->get_context_data(x.get_parms_id())->small_ntt_tables());
    return x;
}

template<typename E, typename R, typename A>
E Verifier<E, R, A>::encode(void *pk, A x) {
    SEALContext *context = static_cast<SEALContext *>(pk);
    Plaintext plain(seal::util::uint_to_hex_string(&x, 1));
    E res = E(*context, plain, &(context->first_parms_id()));
    res.ntt_inplace(context->get_context_data(res.get_parms_id())->small_ntt_tables());
    return res;
}

template<typename E, typename R, typename A>
R Verifier<E, R, A>::decode(void *sk, E x) {
    return x;
}