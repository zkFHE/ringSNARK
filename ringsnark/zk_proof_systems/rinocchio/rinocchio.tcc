#include <ringsnark/reductions/r1cs_to_qrp/r1cs_to_qrp.hpp>

namespace ringsnark::rinocchio {
    template<typename RingT, typename EncT>
    keypair<RingT, EncT> generator(const r1cs_constraint_system<RingT> &cs) {
        const RingT s = RingT::random_exceptional_element();
        const qrp_instance_evaluation<RingT> qrp_inst = r1cs_to_qrp_instance_map_with_evaluation(cs, s);

        const auto [pk_enc, sk_enc] = EncT::keygen();

        const RingT alpha = RingT::random_invertible_element(),
                r_v = RingT::random_invertible_element(),
                r_w = RingT::random_invertible_element(),
                r_y = r_v * r_w,
                beta = RingT::random_nonzero_element();

        // ({E(s^i)}_{i=0}^{num_mid}}, {E(alpha * s^i)}_{i=0}^{num_mid}}, {beta_prod}_{i=0}^{num_mid}, pk)
        // Ht holds the monomials {s^i}_{i=0}^m

        vector<RingT> s_pows_ring(qrp_inst.Ht.begin(), qrp_inst.Ht.begin() + cs.num_constraints() + 1);
        vector<RingT> alpha_s_pows_ring(s_pows_ring);
        for (auto &s_i: alpha_s_pows_ring) { s_i *= alpha; }

        vector<RingT> linchecks, rv_vs, rw_ws, ry_ys;
        linchecks.reserve(cs.auxiliary_input_size);
        rv_vs.reserve(cs.auxiliary_input_size);
        rw_ws.reserve(cs.auxiliary_input_size);
        ry_ys.reserve(cs.auxiliary_input_size);
        for (size_t i = 0; i < cs.auxiliary_input_size; i++) {
            size_t idx = i + cs.primary_input_size + 1;
            rv_vs.push_back(r_v * qrp_inst.At[idx]);
            rw_ws.push_back(r_w * qrp_inst.Bt[idx]);
            ry_ys.push_back(r_y * qrp_inst.Ct[idx]);

            RingT lincheck(rv_vs[i]);
            lincheck += rw_ws[i];
            lincheck += ry_ys[i];
            lincheck *= beta;
            linchecks.push_back(lincheck);
        }

        const vector<EncT> s_pows = EncT::encode(sk_enc, s_pows_ring),
                alpha_s_pows = EncT::encode(sk_enc, alpha_s_pows_ring),
                beta_prods = EncT::encode(sk_enc, linchecks);

        const RingT beta_Zt = beta * qrp_inst.Zt;
        const EncT beta_rv_ts = EncT::encode(sk_enc, {beta_Zt * r_v})[0];
        const EncT beta_rw_ts = EncT::encode(sk_enc, {beta_Zt * r_w})[0];
        const EncT beta_ry_ts = EncT::encode(sk_enc, {beta_Zt * r_y})[0];

        const RingT alpha_Zt = alpha * qrp_inst.Zt;
        const EncT alpha_rv_ts = EncT::encode(sk_enc, {alpha_Zt * r_v})[0];
        const EncT alpha_rw_ts = EncT::encode(sk_enc, {alpha_Zt * r_w})[0];
        const EncT alpha_ry_ts = EncT::encode(sk_enc, {alpha_Zt * r_y})[0];

        // pk = ({E(s^i)}_{i=0}^{num_mid}}, {E(alpha * s^i)}_{i=0}^{num_mid}}, {beta_prod}_{i=0}^{num_mid}, pk_enc)
        auto pk = new proving_key<RingT, EncT>(cs, s_pows, alpha_s_pows, beta_prods,
                                               beta_rv_ts, beta_rw_ts, beta_ry_ts,
                                               alpha_rv_ts, alpha_rw_ts, alpha_ry_ts,
                                               EncT::encode(sk_enc, rv_vs),
                                               EncT::encode(sk_enc, rw_ws),
                                               EncT::encode(sk_enc, ry_ys),
                                               pk_enc);

        // vk = (pk, s, alpha, beta, r_v, r_w, r_y, sk_enc)
        auto vk = new verification_key<RingT, EncT>(*pk, s, alpha, beta, r_v, r_w, r_y, sk_enc);

        return ::ringsnark::rinocchio::keypair<RingT, EncT>(*pk, *vk);
    }

    template<typename RingT, typename EncT>
    proof<RingT, EncT> prover(const proving_key<RingT, EncT> &pk,
                              const r1cs_primary_input<RingT> &primary_input,
                              const r1cs_auxiliary_input<RingT> &auxiliary_input) {
#ifdef DEBUG
        assert(pk.constraint_system.is_satisfied(primary_input, auxiliary_input));
#endif
        const bool use_zk = !auxiliary_input.empty();
        if (!use_zk) {
            cout << "[Prover] " << "using non-zero-knowledge SNARK, since no auxiliary inputs are defined" << endl;
        }
        const RingT d1 = use_zk ? RingT::random_invertible_element() : RingT::zero();
        const RingT d2 = use_zk ? RingT::random_invertible_element() : RingT::zero();
        const RingT d3 = use_zk ? RingT::random_invertible_element() : RingT::zero();


        const qrp_witness<RingT> qrp_wit = r1cs_to_qrp_witness_map(pk.constraint_system,
                                                                   primary_input, auxiliary_input,
                                                                   d1, d2, d3);

        // TODO: this is highly non-optimized, skip all the zero-multiplication
        // s_pows, alpha_s_pows have length d+1, where d = cs.num_constraints() is the size of the QRP
        const auto a_mid = qrp_wit.coefficients_for_A_mid;
        EncT a_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                                a_mid.begin(), a_mid.end());
        EncT alpha_a_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end() - 1,
                                                      a_mid.begin(), a_mid.end());

        const auto b_mid = qrp_wit.coefficients_for_B_mid;
        EncT b_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                                b_mid.begin(), b_mid.end());
        EncT alpha_b_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end() - 1,
                                                      b_mid.begin(), b_mid.end());

        const auto c_mid = qrp_wit.coefficients_for_C_mid;
        EncT c_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                                c_mid.begin(), c_mid.end());
        EncT alpha_c_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end() - 1,
                                                      c_mid.begin(), c_mid.end());

        const auto z = qrp_wit.coefficients_for_Z;
        EncT z_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(),
                                                z.begin(), z.end());
        EncT alpha_z_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(),
                                                      z.begin(), z.end());

        // Add shift terms
        // TODO: add terms to coefficients_for_{A, B, C} directly, similarly to H
        a_enc += d1 * z_enc;
        alpha_a_enc += d1 * alpha_z_enc;
        b_enc += d2 * z_enc;
        alpha_b_enc += d2 * alpha_z_enc;
        c_enc += d3 * z_enc;
        alpha_c_enc += d3 * alpha_z_enc;

        const auto h = qrp_wit.coefficients_for_H;
        EncT d_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end(),
                                                h.begin(), h.end());
        EncT alpha_d_enc = inner_product<EncT, RingT>(pk.alpha_s_pows.begin(), pk.alpha_s_pows.end(),
                                                      h.begin(), h.end());

        EncT f_enc = inner_product<EncT, RingT>(pk.beta_prods.begin(), pk.beta_prods.end(),
                                                auxiliary_input.begin(), auxiliary_input.end());
        f_enc += d1 * pk.beta_rv_ts;
        f_enc += d2 * pk.beta_rw_ts;
        f_enc += d3 * pk.beta_ry_ts;

        return ringsnark::rinocchio::proof<RingT, EncT>(a_enc, alpha_a_enc,
                                                        b_enc, alpha_b_enc,
                                                        c_enc, alpha_c_enc,
                                                        d_enc, alpha_d_enc,
                                                        f_enc);
    }

    template<typename RingT, typename EncT>
    bool verifier(const verification_key<RingT, EncT> &vk,
                  const r1cs_primary_input<RingT> &primary_input,
                  const proof<RingT, EncT> &proof) {
        const RingT V_mid = EncT::decode(vk.sk_enc, proof.A),
                V_mid_prime = EncT::decode(vk.sk_enc, proof.A_prime),
                W_mid = EncT::decode(vk.sk_enc, proof.B),
                W_mid_prime = EncT::decode(vk.sk_enc, proof.B_prime),
                Y_mid = EncT::decode(vk.sk_enc, proof.C),
                Y_mid_prime = EncT::decode(vk.sk_enc, proof.C_prime),
                H = EncT::decode(vk.sk_enc, proof.D),
                H_prime = EncT::decode(vk.sk_enc, proof.D_prime),
                L_beta = EncT::decode(vk.sk_enc, proof.F);
        // TODO: define one version that is amenable to parallelization (close to current implementation), and an "online" version that re-uses the same object in a decode-check loop

        const auto cs = vk.pk.constraint_system;

        qrp_instance_evaluation<RingT> qrp_inst_eval = r1cs_to_qrp_instance_map_with_evaluation(cs, vk.s);

        // L = beta * ((vk.r_v * V_mid) + (vk.r_w * W_mid) + (vk.r_y * Y_mid))
        RingT L = V_mid * vk.r_v;
        L += W_mid * vk.r_w;
        L += Y_mid * vk.r_y;
        L *= vk.beta;

        // TODO: make this more efficient, and skip all zero-mults
        vector<RingT> padded_primary_assignment(primary_input);
        vector<RingT> zeros(vk.pk.constraint_system.auxiliary_input_size, RingT::zero());
        padded_primary_assignment.insert(padded_primary_assignment.end(), zeros.begin(), zeros.end());
        // TODO: or use the {At, Bt, Ct} members from qrp_inst_eval?
        vector<RingT> v_io(cs.num_constraints()), w_io(cs.num_constraints()), y_io(cs.num_constraints());
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            v_io[i] = cs.constraints[i].a.evaluate(padded_primary_assignment);
            w_io[i] = cs.constraints[i].b.evaluate(padded_primary_assignment);
            y_io[i] = cs.constraints[i].c.evaluate(padded_primary_assignment);
        }
        auto domain = get_evaluation_domain<RingT>(cs.num_constraints());
        vector<RingT> xs(domain->m);
        for (size_t i = 0; i < domain->m; i++) { xs[i] = domain->get_domain_element(i); }
        v_io = interpolate(xs, v_io);
        w_io = interpolate(xs, w_io);
        y_io = interpolate(xs, y_io);

        auto v_io_s = eval(v_io, vk.s);
        auto w_io_s = eval(w_io, vk.s);
        auto y_io_s = eval(y_io, vk.s);

        // P = (v_io(s) + V_mid) * (w_io(s) + W_mid) - (y_io(s) + Y_mid)
        RingT P = V_mid + v_io_s;
        P *= W_mid + w_io_s;
        P -= Y_mid + y_io_s;

        RingT tmp;
        // CHECK: V'_mid = alpha * V_mid
        tmp = V_mid * vk.alpha;
        bool res = true;
        if (V_mid_prime != tmp) {
            res = false;
        }
        // CHECK: W'_mid = alpha * W_mid
        tmp = W_mid * vk.alpha;
        if (W_mid_prime != tmp) {
            res = false;
        }
        // CHECK: Y'_mid = alpha * Y_mid
        tmp = Y_mid * vk.alpha;
        if (Y_mid_prime != tmp) {
            res = false;
        }
        // CHECK: H' = alpha * H
        tmp = H * vk.alpha;
        if (H_prime != tmp) {
            res = false;
        }
        // CHECK: L_beta = L
        if (L != L_beta) {
            res = false;
        }

        tmp = H * qrp_inst_eval.Zt;
        if (P != tmp) {
            res = false;
        }
        return res;
    }
}