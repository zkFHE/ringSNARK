/** @file
 *****************************************************************************
 Declaration of interfaces for evaluation domains.
 Roughly, given a desired size m for the domain, the constructor selects
 a choice of domain S with size ~m that has been selected so to optimize
 - computations of Lagrange polynomials, and
 - FFT/iFFT computations.
 An evaluation domain also provides other other functions, e.g., accessing
 individual elements in S or evaluating its vanishing polynomial.
 The descriptions below make use of the definition of a *Lagrange polynomial*,
 which we recall. Given a field F, a subset S=(a_i)_i of F, and an index idx
 in {0,...,|S-1|}, the idx-th Lagrange polynomial (wrt to subset S) is defined to be
 \f[   L_{idx,S}(z) := prod_{k \neq idx} (z - a_k) / prod_{k \neq idx} (a_{idx} - a_k)   \f]
 Note that, by construction:
 \f[   \forall j \neq idx: L_{idx,S}(a_{idx}) = 1  \text{ and }  L_{idx,S}(a_j) = 0   \f]
 *****************************************************************************
 * @author     This file is part of libfqfft, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EVALUATION_DOMAIN_HPP_
#define EVALUATION_DOMAIN_HPP_

#include <vector>
#include <memory>

using namespace std;

namespace ringsnark {

/**
 * An evaluation domain.
 */
    template<typename RingT>
    class evaluation_domain {
    protected:
        vector<RingT> values;

    public:
        const size_t m;

        /**
         * Construct an evaluation domain S of size m, if possible.
         *
         * (See the function get_evaluation_domain below.)
         */
        explicit evaluation_domain(size_t m);

        /**
         * Get the idx-th element in S.
         */
        RingT get_domain_element(const size_t idx) const;

//        /**
//         * Compute the FFT, over the domain S, of the vector a.
//         */
//        virtual void FFT(std::vector<RingT> &a) = 0;
//
//        /**
//         * Compute the inverse FFT, over the domain S, of the vector a.
//         */
//        virtual void iFFT(std::vector<RingT> &a) = 0;
//
//        /**
//         * Compute the FFT, over the domain g*S, of the vector a.
//         */
//        virtual void cosetFFT(std::vector<RingT> &a, const RingT &g) = 0;
//
//        /**
//         * Compute the inverse FFT, over the domain g*S, of the vector a.
//         */
//        virtual void icosetFFT(std::vector<RingT> &a, const RingT &g) = 0;

        /**
         * Evaluate all Lagrange polynomials.
         *
         * The inputs are:
         * - an integer m
         * - an element t
         * The output is a vector (b_{0},...,b_{m-1})
         * where b_{i} is the evaluation of L_{i,S}(z) at z = t.
         */
        std::vector<RingT> evaluate_all_lagrange_polynomials(const RingT &t) const;

        /**
         * Evaluate the vanishing polynomial of S at the ring element t.
         */
        RingT compute_vanishing_polynomial(const RingT &t) const;

        /**
        * Compute the coefficients of the vanishing polynomial of S.
        */
        vector<RingT> vanishing_polynomial() const;

        /**
         * Add the coefficients of the vanishing polynomial of S to the coefficients of the polynomial H.
         */
        void add_poly_Z(const RingT &coeff, std::vector<RingT> &H) const;

        /**
         * Multiply by the evaluation, on a coset of S, of the inverse of the vanishing polynomial of S.
         */
        void divide_by_Z_on_coset(std::vector<RingT> &P) const;
    };


    template<typename RingT>
    std::shared_ptr<evaluation_domain<RingT> > get_evaluation_domain(size_t min_size);

    int test = 0;
}

#include "evaluation_domain.tcc"

#endif // EVALUATION_DOMAIN_HPP_
