/** @file
 *****************************************************************************

 Implementation of interfaces for a R1CS-to-QRP reduction.

 See r1cs_to_qrp.hpp .

 *****************************************************************************/

#ifndef R1CS_TO_QRP_TCC_
#define R1CS_TO_QRP_TCC_

#include "evaluation_domain.hpp"

namespace ringsnark {

/**
 * Instance map for the R1CS-to-QRP reduction.
 *
 * Namely, given a R1CS constraint system cs, construct a QRP instance for which:
 *   A := (A_0(z),A_1(z),...,A_m(z))
 *   B := (B_0(z),B_1(z),...,B_m(z))
 *   C := (C_0(z),C_1(z),...,C_m(z))
 * where
 *   m = number of variables of the QRP
 * and
 *   each A_i,B_i,C_i is expressed in the Lagrange basis.
 */
    template<typename RingT>
    qrp_instance<RingT> r1cs_to_qrp_instance_map(const r1cs_constraint_system<RingT> &cs) {
        const std::shared_ptr<evaluation_domain<RingT>> domain = get_evaluation_domain<RingT>(
                cs.num_constraints() + cs.num_inputs() + 1);

        std::vector<std::map<size_t, RingT>> A_in_Lagrange_basis(cs.num_variables() + 1);
        std::vector<std::map<size_t, RingT>> B_in_Lagrange_basis(cs.num_variables() + 1);
        std::vector<std::map<size_t, RingT>> C_in_Lagrange_basis(cs.num_variables() + 1);

        /**
         * add and process the constraints
         *     input_i * 0 = 0
         * to ensure soundness of input consistency
         */
        // TODO: check if this is also relevant in the ring, if this should be removed, or if it should be handled at a different abstraction level (e.g., in the R1CS generation)
        for (size_t i = 0; i <= cs.num_inputs(); ++i) {
            A_in_Lagrange_basis[i][cs.num_constraints() + i] = RingT::one();
        }
        /* process all other constraints */
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j) {
                A_in_Lagrange_basis[cs.constraints[i].a.terms[j].index][i] += cs.constraints[i].a.terms[j].coeff;
            }

            for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j) {
                B_in_Lagrange_basis[cs.constraints[i].b.terms[j].index][i] += cs.constraints[i].b.terms[j].coeff;
            }

            for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j) {
                C_in_Lagrange_basis[cs.constraints[i].c.terms[j].index][i] += cs.constraints[i].c.terms[j].coeff;
            }
        }

        return qrp_instance<RingT>(domain,
                                   cs.num_variables(),
                                   domain->m,
                                   cs.num_inputs(),
                                   std::move(A_in_Lagrange_basis),
                                   std::move(B_in_Lagrange_basis),
                                   std::move(C_in_Lagrange_basis));
    }


/**
 * Instance map for the R1CS-to-QRP reduction followed by evaluation of the resulting QRP instance.
 *
 * Namely, given a R1CS constraint system cs and a field element t, construct
 * a QRP instance (evaluated at t) for which:
 *   At := (A_0(t),A_1(t),...,A_m(t))
 *   Bt := (B_0(t),B_1(t),...,B_m(t))
 *   Ct := (C_0(t),C_1(t),...,C_m(t))
 *   Ht := (1,t,t^2,...,t^n)
 *   Zt := Z(t) = "vanishing polynomial of a certain set S, evaluated at t"
 * where
 *   m = number of variables of the QRP
 *   n = degree of the QRP
 */
    template<typename RingT>
    qrp_instance_evaluation<RingT> r1cs_to_qrp_instance_map_with_evaluation(const r1cs_constraint_system<RingT> &cs,
                                                                            const RingT &t) {
        const std::shared_ptr<evaluation_domain<RingT>> domain = get_evaluation_domain<RingT>(
                cs.num_constraints() + cs.num_inputs() + 1
        );

        std::vector<RingT> At, Bt, Ct, Ht;

        At.resize(cs.num_variables() + 1, RingT::zero());
        Bt.resize(cs.num_variables() + 1, RingT::zero());
        Ct.resize(cs.num_variables() + 1, RingT::zero());
        Ht.reserve(domain->m + 1);

        const RingT Zt = domain->compute_vanishing_polynomial(t);

        const std::vector<RingT> u = domain->evaluate_all_lagrange_polynomials(t);
        /**
         * add and process the constraints
         *     input_i * 0 = 0
         * to ensure soundness of input consistency
         */
        for (size_t i = 0; i <= cs.num_inputs(); ++i) {
            At[i] = u[cs.num_constraints() + i];
        }
        /* process all other constraints */
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            for (size_t j = 0; j < cs.constraints[i].a.terms.size(); ++j) {
                At[cs.constraints[i].a.terms[j].index] += u[i] * cs.constraints[i].a.terms[j].coeff;
            }

            for (size_t j = 0; j < cs.constraints[i].b.terms.size(); ++j) {
                Bt[cs.constraints[i].b.terms[j].index] += u[i] * cs.constraints[i].b.terms[j].coeff;
            }

            for (size_t j = 0; j < cs.constraints[i].c.terms.size(); ++j) {
                Ct[cs.constraints[i].c.terms[j].index] += u[i] * cs.constraints[i].c.terms[j].coeff;
            }
        }

        RingT ti = RingT::one();
        for (size_t i = 0; i < domain->m + 1; ++i) {
            Ht.emplace_back(ti);
            ti *= t;
        }
        // libff::leave_block("Compute evaluations of A, B, C, H at t");

        // libff::leave_block("Call to r1cs_to_qrp_instance_map_with_evaluation");

        return qrp_instance_evaluation<RingT>(domain,
                                              cs.num_variables(),
                                              domain->m,
                                              cs.num_inputs(),
                                              t,
                                              std::move(At),
                                              std::move(Bt),
                                              std::move(Ct),
                                              std::move(Ht),
                                              Zt);
    }



    // TODO: update description
/**
 * Witness map for the R1CS-to-QRP reduction.
 *
 * The witness map takes zero knowledge into account when d1,d2,d3 are random.
 *
 * More precisely, compute the coefficients
 *     h_0,h_1,...,h_n
 * of the polynomial
 *     H(z) := (A(z)*B(z)-C(z))/Z(z)
 * where
 *   A(z) := A_0(z) + \sum_{k=1}^{m} w_k A_k(z) + d1 * Z(z)
 *   B(z) := B_0(z) + \sum_{k=1}^{m} w_k B_k(z) + d2 * Z(z)
 *   C(z) := C_0(z) + \sum_{k=1}^{m} w_k C_k(z) + d3 * Z(z)
 *   Z(z) := "vanishing polynomial of set S"
 * and
 *   m = number of variables of the QRP
 *   n = degree of the QRP
 *
 * This is done as follows:
 *  (1) compute evaluations of A,B,C on S = {sigma_1,...,sigma_n}
 *  (2) compute coefficients of A,B,C
 *  (3) compute evaluations of A,B,C on T = "coset of S"
 *  (4) compute evaluation of H on T
 *  (5) compute coefficients of H
 *  (6) patch H to account for d1,d2,d3 (i.e., add coefficients of the polynomial (A d2 + B d1 - d3) + d1*d2*Z )
 *
 * The code below is not as simple as the above high-level description due to
 * some reshuffling to save space.
 */
    template<typename RingT>
    qrp_witness<RingT> r1cs_to_qrp_witness_map(const r1cs_constraint_system<RingT> &cs,
                                               const r1cs_primary_input<RingT> &primary_input,
                                               const r1cs_auxiliary_input<RingT> &auxiliary_input,
                                               const RingT &d1,
                                               const RingT &d2,
                                               const RingT &d3) {
        /* sanity check */
        assert(cs.is_satisfied(primary_input, auxiliary_input));

        const std::shared_ptr<evaluation_domain<RingT>>
                domain = get_evaluation_domain<RingT>(
                cs.num_constraints() + cs.num_inputs() + 1);

        r1cs_variable_assignment<RingT> full_variable_assignment = primary_input;
        full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());


        // Compute coefficiens for A_mid, B_mid, C_mid
        std::vector<RingT> a_mid(domain->m, RingT::zero()), b_mid(domain->m, RingT::zero()), c_mid(domain->m,
                                                                                                   RingT::zero());
        for (size_t i = 0; i <= cs.num_inputs(); ++i) {
            a_mid[i + cs.num_constraints()] = (i > 0 ? full_variable_assignment[i - 1] : RingT::one());
        }
        /* account for all other constraints */
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            a_mid[i] += cs.constraints[i].a.evaluate(auxiliary_input);
            b_mid[i] += cs.constraints[i].b.evaluate(auxiliary_input);
            c_mid[i] += cs.constraints[i].c.evaluate(auxiliary_input);
        }

        // Compute coefficients for vanishing polynomial Z
        std::vector<RingT> Z = domain->vanishing_polynomial();


        // Compute coefficients for H
        std::vector<RingT> aA(domain->m, RingT::zero()), aB(domain->m, RingT::zero()), aC(domain->m, RingT::zero());

        for (size_t i = 0; i <= cs.num_inputs(); ++i) {
            aA[i + cs.num_constraints()] = (i > 0 ? full_variable_assignment[i - 1] : RingT::one());
        }
        /* account for all other constraints */
        for (size_t i = 0; i < cs.num_constraints(); ++i) {
            aA[i] += cs.constraints[i].a.evaluate(full_variable_assignment);
            aB[i] += cs.constraints[i].b.evaluate(full_variable_assignment);
            aC[i] += cs.constraints[i].c.evaluate(full_variable_assignment);
        }

//        domain->iFFT(aA);
//
//        domain->iFFT(aB);

        std::vector<RingT> coefficients_for_H(domain->m + 1, RingT::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
        /* add coefficients of the polynomial (d2*A + d1*B - d3) + d1*d2*Z */
        for (size_t i = 0; i < domain->m; ++i) {
            coefficients_for_H[i] = d2 * aA[i] + d1 * aB[i];
        }
        coefficients_for_H[0] -= d3;
        domain->add_poly_Z(d1 * d2, coefficients_for_H);

//        domain->cosetFFT(aA, RingT::multiplicative_generator);
//
//        domain->cosetFFT(aB, RingT::multiplicative_generator);

        // Compute coefficients of A*B - C
        // TODO: can we use 2*m-1 instead (even for FFT impl. as some later point)?
        auto m = cs.num_constraints() + cs.num_inputs() + 1;
        auto domain_2m = get_evaluation_domain<RingT>(2 * m);
        vector<RingT> x, vwy;
        x.reserve(2 * m);
        vwy.reserve(2 * m);
        for (size_t i = 0; i < 2 * m; i++) {
            x.push_back(domain_2m->get_domain_element(i));
            vector<RingT> v(m, RingT::zero()), w(m, RingT::zero()), y(m, RingT::zero());
            RingT vwy_i = eval(aA, x[i]);
            vwy_i *= eval(aB, x[i]);
            vwy_i -= eval(aC, x[i]);

            vwy.push_back(vwy_i);
        }
        vector<RingT> H_tmp = interpolate(x, vwy);
        domain->divide_by_Z_on_coset(H_tmp);
        assert(coefficients_for_H.size() == H_tmp.size());

//        std::vector<RingT> &H_tmp = aA; // can overwrite aA because it is not used later
//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
//        for (size_t i = 0; i < domain->m; ++i) {
//            H_tmp[i] = aA[i] * aB[i];
//        }
//        std::vector<RingT>().swap(aB); // destroy aB


//        domain->iFFT(aC);
//
//        domain->cosetFFT(aC, RingT::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
//        for (size_t i = 0; i < domain->m; ++i) {
//            H_tmp[i] = (H_tmp[i] - aC[i]);
//        }

//        domain->divide_by_Z_on_coset(H_tmp);
//
//        domain->icosetFFT(H_tmp, RingT::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
        for (size_t i = 0; i < domain->m; ++i) {
            coefficients_for_H[i] += H_tmp[i];
        }


        return qrp_witness<RingT>(cs.num_variables(),
                                  domain->m,
                                  cs.num_inputs(),
                                  d1,
                                  d2,
                                  d3,
                                  full_variable_assignment,
                                  a_mid, b_mid, c_mid, Z,
                                  std::move(coefficients_for_H));
    }

} // ringsnark

#endif // R1CS_TO_QRP_TCC_