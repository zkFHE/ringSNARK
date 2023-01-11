/** @file
 *****************************************************************************

 Declaration of interfaces for a QRP ("Quadratic Ring Program").

 QRPs are defined in \[GNSV21], and are very similar to QAPs (Quadratic Arithmetic Programs) \[GGPR13].

 References:

 \[GNSV21]:
 "Rinocchio: SNARKs for Ring Arithmetic",
 Chaya Ganesh, Anca Nitulescu, Eduardo Soria-Vazquez,
<https://eprint.iacr.org/2021/322>

 \[GGPR13]:
 "Quadratic span programs and succinct NIZKs without PCPs",
 Rosario Gennaro, Craig Gentry, Bryan Parno, Mariana Raykova,
 EUROCRYPT 2013,
 <http://eprint.iacr.org/2012/215>

 *****************************************************************************/

#ifndef QRP_HPP_
#define QRP_HPP_

#include <map>
#include <memory>

#include "ringsnark/util/evaluation_domain.hpp"

namespace ringsnark {

/* forward declaration */
    template<typename RingT>
    class qrp_witness;

/**
 * A QRP instance.
 *
 * Specifically, the datastructure stores:
 * - a choice of domain (corresponding to a certain subset of the ring);
 * - the number of variables, the degree, and the number of inputs; and
 * - coefficients of the A,B,C polynomials in the Lagrange basis.
 *
 * There is no need to store the Z polynomial because it is uniquely
 * determined by the domain (as Z is its vanishing polynomial).
 */
    template<typename RingT>
    class qrp_instance {
    private:
        size_t num_variables_;
        size_t degree_;
        size_t num_inputs_;

    public:
        std::shared_ptr<evaluation_domain<RingT> > domain;

        std::vector<std::map<size_t, RingT> > A_in_Lagrange_basis;
        std::vector<std::map<size_t, RingT> > B_in_Lagrange_basis;
        std::vector<std::map<size_t, RingT> > C_in_Lagrange_basis;

        std::vector<std::map<size_t, RingT> > A;
        std::vector<std::map<size_t, RingT> > B;
        std::vector<std::map<size_t, RingT> > C;

        qrp_instance(const std::shared_ptr<evaluation_domain<RingT> > &domain,
                     size_t num_variables,
                     size_t degree,
                     size_t num_inputs,
                     const std::vector<std::map<size_t, RingT> > &A_in_Lagrange_basis,
                     const std::vector<std::map<size_t, RingT> > &B_in_Lagrange_basis,
                     const std::vector<std::map<size_t, RingT> > &C_in_Lagrange_basis
        );

        qrp_instance(const std::shared_ptr<evaluation_domain<RingT> > &domain,
                     size_t num_variables,
                     size_t degree,
                     size_t num_inputs,
                     std::vector<std::map<size_t, RingT> >
                     &&A_in_Lagrange_basis,
                     std::vector<std::map<size_t, RingT> > &&B_in_Lagrange_basis,
                     std::vector<std::map<size_t, RingT> >
                     &&C_in_Lagrange_basis);

        qrp_instance(const qrp_instance<RingT> &other) = default;

        qrp_instance(qrp_instance<RingT> &&other) noexcept = default;

        qrp_instance &operator=(const qrp_instance<RingT> &other) = default;

        qrp_instance &operator=(qrp_instance<RingT> &&other) noexcept = default;

        [[nodiscard]] size_t num_variables() const;

        [[nodiscard]] size_t degree() const;

        [[nodiscard]] size_t num_inputs() const;

        bool is_satisfied(const qrp_witness<RingT> &witness) const;
    };

/**
 * A QRP instance evaluation is a QRP instance that is evaluated at a ring element t.
 *
 * Specifically, the datastructure stores:
 * - a choice of domain (corresponding to a certain subset of the ring);
 * - the number of variables, the degree, and the number of inputs;
 * - a ring element t;
 * - evaluations of the A,B,C (and Z) polynomials at t;
 * - evaluations of all monomials of t;
 * - counts about how many of the above evaluations are in fact non-zero.
 */
    template<typename RingT>
    class qrp_instance_evaluation {
    private:
        size_t num_variables_;
        size_t degree_;
        size_t num_inputs_;
    public:
        std::shared_ptr<evaluation_domain<RingT> > domain;

        RingT t;

        std::vector<RingT> At, Bt, Ct, Ht;
        std::vector<RingT> Aio_t, Bio_t, Cio_t;

        RingT Zt;

        qrp_instance_evaluation(const std::shared_ptr<evaluation_domain<RingT> > &domain,
                                size_t num_variables,
                                size_t degree,
                                size_t num_inputs,
                                const RingT &t,
                                const std::vector<RingT> &At,
                                const std::vector<RingT> &Bt,
                                const std::vector<RingT> &Ct,
                                const std::vector<RingT> &Ht,
                                const RingT &Zt
        );

        qrp_instance_evaluation(const std::shared_ptr<evaluation_domain<RingT> > &domain,
                                size_t num_variables,
                                size_t degree,
                                size_t num_inputs,
                                const RingT &t,
                                std::vector<RingT>
                                &&At,
                                std::vector<RingT> &&Bt,
                                std::vector<RingT>
                                &&Ct,
                                std::vector<RingT> &&Ht,
                                const RingT &Zt
        );

        qrp_instance_evaluation(const qrp_instance_evaluation<RingT> &other) = default;

        qrp_instance_evaluation(qrp_instance_evaluation<RingT> &&other) noexcept = default;

        qrp_instance_evaluation &operator=(const qrp_instance_evaluation<RingT> &other) = default;

        qrp_instance_evaluation &operator=(qrp_instance_evaluation<RingT> &&other) noexcept = default;

        [[nodiscard]] size_t num_variables() const;

        [[nodiscard]] size_t degree() const;

        [[nodiscard]] size_t num_inputs() const;

        bool is_satisfied(const qrp_witness<RingT> &witness) const;
    };

/**
 * A QRP witness.
 */
    template<typename RingT>
    class qrp_witness {
    private:
        const size_t num_variables_;
        const size_t degree_;
        const size_t num_inputs_;

    public:
        const RingT d1, d2, d3;

        const std::vector<RingT> coefficients_for_ABCs;
        const std::vector<RingT> coefficients_for_A_io;
        const std::vector<RingT> coefficients_for_B_io;
        const std::vector<RingT> coefficients_for_C_io;
        const std::vector<RingT> coefficients_for_A_mid;
        const std::vector<RingT> coefficients_for_B_mid;
        const std::vector<RingT> coefficients_for_C_mid;
        const std::vector<RingT> coefficients_for_Z;
        const std::vector<RingT> coefficients_for_H;

        qrp_witness(size_t num_variables,
                    size_t degree,
                    size_t num_inputs,
                    const RingT &d1,
                    const RingT &d2,
                    const RingT &d3,
                    const std::vector<RingT> &coefficients_for_ABCs,
                    const std::vector<RingT> &coefficients_for_A_io,
                    const std::vector<RingT> &coefficients_for_B_io,
                    const std::vector<RingT> &coefficients_for_C_io,
                    const std::vector<RingT> &coefficients_for_A_mid,
                    const std::vector<RingT> &coefficients_for_B_mid,
                    const std::vector<RingT> &coefficients_for_C_mid,
                    const std::vector<RingT> &coefficients_for_Z,
                    const std::vector<RingT> &coefficients_for_H
        );

        qrp_witness(size_t num_variables,
                    size_t degree,
                    size_t num_inputs,
                    const RingT &d1,
                    const RingT &d2,
                    const RingT &d3,
                    const std::vector<RingT> &coefficients_for_ABCs,
                    const std::vector<RingT> &coefficients_for_A_io,
                    const std::vector<RingT> &coefficients_for_B_io,
                    const std::vector<RingT> &coefficients_for_C_io,
                    const std::vector<RingT> &coefficients_for_A_mid,
                    const std::vector<RingT> &coefficients_for_B_mid,
                    const std::vector<RingT> &coefficients_for_C_mid,
                    std::vector<RingT> &&coefficients_for_H);

        qrp_witness(const qrp_witness<RingT> &other) = default;

        qrp_witness(qrp_witness<RingT> &&other) noexcept = default;

        qrp_witness &operator=(const qrp_witness<RingT> &other) = default;

        qrp_witness &operator=(qrp_witness<RingT> &&other) noexcept = default;

        [[nodiscard]] size_t num_variables() const;

        [[nodiscard]] size_t degree() const;

        [[nodiscard]] size_t num_inputs() const;
    };

} // ringsnark

#include <ringsnark/relations/arithmetic_programs/qrp/qrp.tcc>

#endif // QRP_HPP_
