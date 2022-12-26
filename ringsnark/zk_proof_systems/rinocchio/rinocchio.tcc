#include <ringsnark/reductions/r1cs_to_qrp/r1cs_to_qrp.hpp>

namespace ringsnark::rinocchio {
    template<typename RingT, typename EncT>
    keypair<RingT, EncT> generator(const r1cs_constraint_system<RingT> &cs) {
        RingT s = RingT::random_exceptional_element();

        qrp_instance<RingT> qrp_instance = r1cs_to_qrp_instance_map(cs);
        qrp_instance_evaluation<RingT> qrp_inst = r1cs_to_qrp_instance_map_with_evaluation(cs, s);

        auto [pk_enc, sk_enc] = EncT::keygen();

        const RingT alpha = RingT::random_invertible_element(),
                r_v = RingT::random_invertible_element(),
                r_w = RingT::random_invertible_element(),
                r_y = r_v * r_w,
                beta = RingT::random_nonzero_element();

        // ({E(s^i)}_{i=0}^{num_mid}}, {E(alpha * s^i)}_{i=0}^{num_mid}}, {beta_prod}_{i=0}^{num_mid}, pk)
        // Ht holds the monomials {s^i}_{i=0}^m

        vector<RingT> alpha_s_pows_ring(qrp_inst.Ht);
        for (auto &s_i: alpha_s_pows_ring) {
            s_i *= alpha;
        }

        vector<RingT> linchecks;
        linchecks.reserve(qrp_inst.At.size());
        for (size_t i = 0; i < qrp_inst.At.size(); i++) {
            RingT lincheck = r_v * qrp_inst.At[i];
            lincheck += r_w * qrp_inst.Bt[i];
            lincheck += r_y * qrp_inst.Ct[i];
            lincheck *= beta;
            linchecks.push_back(lincheck);
        }

        vector<EncT> s_pows = EncT::encode(sk_enc, qrp_inst.Ht),
                alpha_s_pows = EncT::encode(sk_enc, alpha_s_pows_ring),
                beta_prods = EncT::encode(sk_enc, linchecks);

        proving_key<RingT, EncT> pk(cs, s_pows, alpha_s_pows, beta_prods, pk_enc);
        verification_key<RingT, EncT> vk(pk, s, alpha, beta, r_v, r_w, r_y, sk_enc);
        assert(pk.constraint_system.constraints.size() == 1);

        ringsnark::rinocchio::keypair<RingT, EncT> keypair(pk, vk);
        assert(keypair.pk.constraint_system.constraints.size() == 1);
        return keypair;
    }

    template<typename RingT, typename EncT>
    proof<RingT, EncT> prover(const proving_key<RingT, EncT> &pk,
                              const r1cs_primary_input<RingT> &primary_input,
                              const r1cs_auxiliary_input<RingT> &auxiliary_input) {
        assert(pk.constraint_system.constraints.size() == 1);


#ifdef DEBUG
        assert(pk.constraint_system.is_satisfied(primary_input, auxiliary_input));
#endif
        // TODO: skip if used only as a SNARK, without ZK?
        const RingT d1 = RingT::random_invertible_element();
        const RingT d2 = RingT::random_invertible_element();
        const RingT d3 = RingT::random_invertible_element();

        const qrp_witness<RingT> qrp_wit = r1cs_to_qrp_witness_map(pk.constraint_system,
                                                                   primary_input, auxiliary_input, d1, d2, d3);

        const vector<RingT> a_mid = qrp_wit.coefficients_for_A_mid,
                b_mid = qrp_wit.coefficients_for_B_mid,
                c_mid = qrp_wit.coefficients_for_C_mid,
                z = qrp_wit.coefficients_for_Z,
                h = qrp_wit.coefficients_for_H;

        EncT a_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(), a_mid.begin(), a_mid.end());
        EncT alpha_a_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(), a_mid.begin(),
                                                      a_mid.end());
        EncT b_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(), b_mid.begin(), b_mid.end());
        EncT alpha_b_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(), b_mid.begin(),
                                                      b_mid.end());
        EncT c_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(), c_mid.begin(), c_mid.end());
        EncT alpha_c_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(), c_mid.begin(),
                                                      c_mid.end());

        EncT z_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.begin() + z.size(), z.begin(), z.end());

        // Add shift terms
        // TODO: add terms to coefficients_for_{A, B, C} directly, similarly to H
        EncT z1(z_enc);
        z1 *= d1;
        EncT z2(z_enc);
        z2 *= d2;
        EncT z3(z_enc);
        z3 *= d3;
        a_enc += z1;
        alpha_a_enc += z1;
        b_enc += z2;
        alpha_b_enc += z2;
        c_enc += z3;
        alpha_c_enc += z3;

        EncT d_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(), h.begin(), h.end());
        EncT alpha_d_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(), h.begin(),
                                                      h.end());

        EncT f_enc(pk.beta_prods[0]);
        for (size_t i = 1; i < pk.beta_prods.size(); i++) {
            f_enc += pk.beta_prods[i];
        }

        return ringsnark::rinocchio::proof<RingT, EncT>(a_enc, alpha_a_enc,
                                                        b_enc, alpha_b_enc,
                                                        c_enc, alpha_c_enc,
                                                        d_enc, alpha_d_enc, f_enc);
    }

    template<typename RingT, typename EncT>
    bool verifier(const verification_key<RingT, EncT> &vk,
                  const r1cs_primary_input<RingT> &primary_input,
                  const proof<RingT, EncT> &proof) {
        RingT V_mid = EncT::decode(vk.sk_enc, proof.A);
        RingT V_mid_prime = EncT::decode(vk.sk_enc, proof.A_prime);
        RingT W_mid = EncT::decode(vk.sk_enc, proof.B);
        RingT W_mid_prime = EncT::decode(vk.sk_enc, proof.B_prime);
        RingT Y_mid = EncT::decode(vk.sk_enc, proof.C);
        RingT Y_mid_prime = EncT::decode(vk.sk_enc, proof.C_prime);
        RingT H = EncT::decode(vk.sk_enc, proof.D);
        RingT H_prime = EncT::decode(vk.sk_enc, proof.D_prime);
        RingT L_beta = EncT::decode(vk.sk_enc, proof.F);



        // A = E(V_mid)
        // A' = E(V'_mid)
        // B = E(W_mid)
        // B' = E(W'_mid)
        // C = E(Y_mid)
        // C' = E(Y'_mid)
        // D = E(H)
        // D' = E(H')
        // F = E(L)


        // v_io(x) = sum_{k \in I_io} a_k * v_k(x)
        // w_io(x) = sum_{k \in I_io} a_k * w_k(x)
        // y_io(x) = sum_{k \in I_io} a_k * y_k(x)
        auto cs = vk.pk.constraint_system;

        qrp_instance_evaluation<RingT> qrp_inst_eval = r1cs_to_qrp_instance_map_with_evaluation(cs, vk.s);

        // L = beta * ((vk.r_v * V_mid) + (vk.r_w * W_mid) + (vk.r_y * Y_mid))
        RingT L = V_mid * vk.r_v;
        L += W_mid * vk.r_w;
        L += Y_mid * vk.r_y;
        L *= vk.beta;

        RingT v_io_s = qrp_inst_eval.At[0], w_io_s = qrp_inst_eval.Bt[0], y_io_s = qrp_inst_eval.Ct[0];
        /* account for all other constraints */
        // TODO: or use the {At, Bt, Ct} members from qrp_inst_eval?
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            v_io_s += cs.constraints[i].a.evaluate(primary_input);
            w_io_s += cs.constraints[i].b.evaluate(primary_input);
            y_io_s += cs.constraints[i].c.evaluate(primary_input);
        }

        // P = (v_io(s) + V_mid) * (w_io(s) + W_mid) - (y_io(s) + Y_mid)
        RingT P = V_mid * v_io_s;
        P *= W_mid + w_io_s;
        P -= Y_mid + y_io_s;

        RingT tmp;
        // CHECK: V'_mid = alpha * V_mid
        tmp = V_mid * vk.alpha;
        if (V_mid_prime != tmp) {
            return false;
        }
        // CHECK: W'_mid = alpha * W_mid
        tmp = W_mid * vk.alpha;
        if (W_mid_prime != tmp) {
            return false;
        }
        // CHECK: Y'_mid = alpha * Y_mid
        tmp = Y_mid * vk.alpha;
        if (Y_mid_prime != tmp) {
            return false;
        }
        // CHECK: H' = alpha * H
        tmp = H * vk.alpha;
        if (H_prime != tmp) {
            return false;
        }
        // CHECK: L_beta = L
        if (L != L_beta) {
            return false;
        }
        // P = H * t(s)
        tmp = H * qrp_inst_eval.Zt;
        if (P != tmp) {
            return false;
        }
        return true;
    }
}