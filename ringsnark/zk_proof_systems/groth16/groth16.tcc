#include <ringsnark/reductions/r1cs_to_qrp/r1cs_to_qrp.hpp>

namespace ringsnark::groth16 {
    template<typename RingT, typename EncT>
    keypair<RingT, EncT> generator(const r1cs_constraint_system<RingT> &cs) {
        const RingT s = RingT::random_exceptional_element();
        const qrp_instance_evaluation<RingT> qrp_inst = r1cs_to_qrp_instance_map_with_evaluation(cs, s);

        const auto [pk_enc, sk_enc] = EncT::keygen();

        const RingT alpha = RingT::random_invertible_element(),
                beta = RingT::random_invertible_element(),
                gamma = RingT::random_invertible_element(),
                delta = RingT::random_invertible_element();
        // TODO: get inverse directly from generation, since it has to be computed anyway?
        const RingT gamma_inv = RingT::one() / gamma;
        const RingT delta_inv = RingT::one() / delta;

        // ({E(s^i)}_{i=0}^{num_mid}}, {E(alpha * s^i)}_{i=0}^{num_mid}}, {beta_prod}_{i=0}^{num_mid}, pk)
        // Ht holds the monomials {s^i}_{i=0}^m

        vector<RingT> s_pows_ring(qrp_inst.Ht.begin(), qrp_inst.Ht.begin() + cs.num_constraints() + 1);
        vector<RingT> delta_ts_ring(s_pows_ring);
        for (auto &s_i: delta_ts_ring) {
            s_i *= qrp_inst.Zt;
            s_i *= delta_inv;
        }

        RingT tmp;

        vector<RingT> gamma_io;
        gamma_io.reserve(cs.primary_input_size + 1);
        for (size_t i = 0; i < cs.primary_input_size + 1; i++) {
            tmp = beta * qrp_inst.At[i];
            tmp += alpha * qrp_inst.Bt[i];
            tmp += qrp_inst.Ct[i];
            tmp *= gamma_inv;
            gamma_io.push_back(tmp);
        }

        vector<RingT> delta_mid;
        delta_mid.reserve(cs.auxiliary_input_size);
        for (size_t i = 0; i < cs.auxiliary_input_size; i++) {
            size_t idx = i + cs.primary_input_size + 1;

            tmp = beta * qrp_inst.At[idx];
            tmp += alpha * qrp_inst.Bt[idx];
            tmp += qrp_inst.Ct[idx];
            tmp *= delta_inv;
            delta_mid.push_back(tmp);
        }

        auto pk = new proving_key<RingT, EncT>(cs,
                                               EncT::encode(sk_enc, {alpha})[0],
                                               EncT::encode(sk_enc, {beta})[0],
                                               EncT::encode(sk_enc, s_pows_ring),
                                               EncT::encode(sk_enc, gamma_io),
                                               EncT::encode(sk_enc, delta_mid),
                                               EncT::encode(sk_enc, delta_ts_ring),
                                               pk_enc);

        auto vk = new verification_key<RingT, EncT>(*pk, s, alpha, beta, gamma, delta, sk_enc);

        return ::ringsnark::groth16::keypair<RingT, EncT>(*pk, *vk);
    }

    template<typename RingT, typename EncT>
    proof<RingT, EncT> prover(const proving_key<RingT, EncT> &pk,
                              const r1cs_primary_input<RingT> &primary_input,
                              const r1cs_auxiliary_input<RingT> &auxiliary_input) {
#ifdef DEBUG
        assert(pk.constraint_system.is_satisfied(primary_input, auxiliary_input));
#endif
        const bool use_zk = false;
        if (!use_zk) {
            cout << "[Prover] " << "using non-zero-knowledge SNARK" << endl;
        }


        const qrp_witness<RingT> qrp_wit = r1cs_to_qrp_witness_map(pk.constraint_system,
                                                                   primary_input, auxiliary_input,
                                                                   RingT::zero(), RingT::zero(), RingT::zero());

        // TODO: this is highly non-optimized, skip all the zero-multiplication
        // s_pows have length d+1, where d = cs.num_constraints() is the size of the QRP
        auto a = qrp_wit.coefficients_for_A_io;
        EncT a_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                                a.begin(), a.end());
        a = qrp_wit.coefficients_for_A_mid;
        a_enc += inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                            a.begin(), a.end());
        a_enc += pk.alpha;

        auto b = qrp_wit.coefficients_for_B_io;
        EncT b_enc = inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                                b.begin(), b.end());
        b = qrp_wit.coefficients_for_B_mid;
        b_enc += inner_product<EncT, RingT>(pk.s_pows.begin(), pk.s_pows.end() - 1,
                                            b.begin(), b.end());
        b_enc += pk.beta;

        auto h = qrp_wit.coefficients_for_H;
        EncT c_enc = inner_product<EncT, RingT>(pk.delta_ts.begin(), pk.delta_ts.end(),
                                                h.begin(), h.end());
        c_enc += inner_product<EncT, RingT>(pk.delta_mid.begin(), pk.delta_mid.end(),
                                            auxiliary_input.begin(), auxiliary_input.end());

        return ringsnark::groth16::proof<RingT, EncT>(a_enc, b_enc, c_enc);
    }

    template<typename RingT, typename EncT>
    bool verifier(const verification_key<RingT, EncT> &vk,
                  const r1cs_primary_input<RingT> &primary_input,
                  const proof <RingT, EncT> &proof) {
        const RingT A = EncT::decode(vk.sk_enc, proof.A),
                B = EncT::decode(vk.sk_enc, proof.B),
                C = EncT::decode(vk.sk_enc, proof.C);

        const auto cs = vk.pk.constraint_system;

        qrp_instance_evaluation<RingT> qrp_inst_eval = r1cs_to_qrp_instance_map_with_evaluation(cs, vk.s);


        // TODO: make this more efficient, and skip all zero-mults
        vector<RingT> padded_primary_assignment(primary_input);
        vector<RingT> zeros(vk.pk.constraint_system.auxiliary_input_size, RingT::zero());
        padded_primary_assignment.insert(padded_primary_assignment.end(), zeros.begin(), zeros.end());
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

        // TODO: do we need to chek on encodings, or can we only do checks on plaintexts/ring elements?
        // If the latter, can  we simplify the gamma * (...) / gamma term in rhs?
        RingT f_io = vk.beta * v_io_s;
        f_io += vk.alpha * w_io_s;
        f_io += y_io_s;
        f_io /= vk.gamma;

        RingT AB = A * B;
        RingT rhs = vk.alpha * vk.beta;
        rhs += vk.gamma * f_io;
        rhs += vk.delta * C;

        return AB == rhs;
    }
}