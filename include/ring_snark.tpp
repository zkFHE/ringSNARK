#include <iostream>
#include <cmath>
#include "ring_snark.h"


using namespace std;
using namespace rinocchio;

namespace rinocchio {

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
        for (size_t k = params.num_io; k < params.num_io + params.num_mid; k++) {
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
        E f = sum(crs.beta_prod);

        return Proof<R>{a, a_alpha, b, b_alpha, c, c_alpha, d, d_alpha, f};
    }

    template<typename E, typename R, typename A>
    E Prover<E, R, A>::sum(const vector<E> &values) {
        E res(values[0]);
        for (size_t j = 1; j < values.size(); j++) {
            res.add_inplace(values[j]);
        }
        return res;
    }

    template<typename E, typename R, typename A>
    E Prover<E, R, A>::inner_prod_enc_ring(const vector<R> &poly, const vector<E> &values) {
        assert(poly.size() <= values.size()); // poly might have lower degree
        E res(values[0]);
        res.multiply_inplace(poly[0]);
        for (size_t j = 1; j < poly.size(); j++) {
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
            E tmp(values[j]);
            tmp.multiply_scalar_inplace(poly[j]);
            res.add_inplace(tmp);
        }
        return res;
    }

    template<typename R, typename A>
    SnarkParameters<R, A> setup(uint64_t security_param) {
        return SnarkParameters<R, A>();
    };

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
        R v_io_s = params.v0.eval_at(vk.s);
        R w_io_s = params.w0.eval_at(vk.s);
        R y_io_s = params.y0.eval_at(vk.s);


        size_t k = 1; // k=0 is handled separately
        for (auto const &c_k: params.values_io) {
            R tmp = R(c_k);
            tmp.multiply_scalar_inplace(eval(params.v[k], vk.s));
            v_io_s.add_inplace(tmp);

            tmp = R(c_k);
            tmp.multiply_scalar_inplace(eval(params.w[k], vk.s));
            w_io_s.add_inplace(tmp);

            tmp = R(c_k);
            tmp.multiply_scalar_inplace(eval(params.y[k], vk.s));
            y_io_s.add_inplace(tmp);

            k++;
        }

        // L_span = r_v * V_mid + r_w * V_mid + r_y * Y_mid
        //R L_span = (vk.r_v * V_mid) + (vk.r_w * W_mid) + (vk.r_y * Y_mid);
        R L = Y_mid;
        L.multiply_inplace(vk.r_y);

        R tmp = W_mid;
        tmp.multiply_inplace(vk.r_w);
        L.add_inplace(tmp);

        tmp = V_mid;
        tmp.multiply_inplace(vk.r_v);
        L.add_inplace(tmp);

        // P = (v_io(s) + V_mid) * (w_io(s) + W_mid) - (y_io(s) + Y_mid)
        // R P = (v_io_s + V_mid) * (w_io_s + W_mid) - (y_io_s + Y_mid);
        R P = V_mid;
        P.add_inplace(v_io_s);

        tmp = W_mid;
        tmp.add_inplace(v_io_s);
        P.multiply_inplace(tmp);

        P.subtract_inplace(y_io_s);
        P.subtract_inplace(Y_mid);

        // CHECK: V'_mid = alpha * V_mid
        tmp = V_mid;
        tmp.multiply_inplace(vk.alpha);
        if (!V_mid_prime.is_equal(tmp)) {
            return false;
        }
        // CHECK: W'_mid = alpha * W_mid
        tmp = W_mid;
        tmp.multiply_inplace(vk.alpha);
        if (!W_mid_prime.is_equal(tmp)) {
            return false;
        }
        // CHECK: Y'_mid = alpha * Y_mid
        tmp = Y_mid;
        tmp.multiply_inplace(vk.alpha);
        if (!Y_mid_prime.is_equal(tmp)) {
            return false;
        }
        // CHECK: H' = alpha * H
        tmp = H;
        tmp.multiply_inplace(vk.alpha);
        if (!H_prime.is_equal(tmp)) {
            return false;
        }
        // CHECK: L = beta * L_span
        tmp = L;
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
        for (int i = poly.size() - 1; i >= 0; i--) {
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
}