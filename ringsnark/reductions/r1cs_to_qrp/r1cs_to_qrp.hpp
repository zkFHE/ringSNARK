/** @file
 *****************************************************************************

 Declaration of interfaces for a R1CS-to-QRP reduction, that is, constructing
 a QRP ("Quadratic Arithmetic Program") from a R1CS ("Rank-1 Constraint System").

 QRPs are defined in \[GGPR13], and constructed for R1CS also in \[GGPR13].

 The implementation of the reduction follows, extends, and optimizes
 the efficient approach described in Appendix E of \[BCGTV13].

 References:

 \[BCGTV13]
 "SNARKs for C: Verifying Program Executions Succinctly and in Zero Knowledge",
 Eli Ben-Sasson, Alessandro Chiesa, Daniel Genkin, Eran Tromer, Madars Virza,
 CRYPTO 2013,
 <http://eprint.iacr.org/2013/507>

 \[GGPR13]:
 "Quadratic span programs and succinct NIZKs without PCPs",
 Rosario Gennaro, Craig Gentry, Bryan Parno, Mariana Raykova,
 EUROCRYPT 2013,
 <http://eprint.iacr.org/2012/215>

 *****************************************************************************/

#ifndef R1CS_TO_QRP_HPP_
#define R1CS_TO_QRP_HPP_

#include <ringsnark/relations/arithmetic_programs/qrp/qrp.hpp>
#include <ringsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>

namespace ringsnark {

/**
 * Instance map for the R1CS-to-QRP reduction.
 */
    template<typename RingT>
    qrp_instance <RingT> r1cs_to_qrp_instance_map(const r1cs_constraint_system<RingT> &cs);

/**
 * Instance map for the R1CS-to-QRP reduction followed by evaluation of the resulting QRP instance.
 */
    template<typename RingT>
    qrp_instance_evaluation <RingT> r1cs_to_qrp_instance_map_with_evaluation(const r1cs_constraint_system<RingT> &cs,
                                                                             const RingT &t);

/**
 * Witness map for the R1CS-to-QRP reduction.
 *
 * The witness map takes zero knowledge into account when d1,d2,d3 are random.
 */
    template<typename RingT>
    qrp_witness <RingT> r1cs_to_qrp_witness_map(const r1cs_constraint_system<RingT> &cs,
                                                const r1cs_primary_input<RingT> &primary_input,
                                                const r1cs_auxiliary_input<RingT> &auxiliary_input,
                                                const RingT &d1,
                                                const RingT &d2,
                                                const RingT &d3);

} // ringsnark

#include <ringsnark/reductions/r1cs_to_qrp/r1cs_to_qrp.tcc>

#endif // R1CS_TO_QRP_HPP_
