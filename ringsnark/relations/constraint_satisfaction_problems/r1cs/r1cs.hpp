/** @file
 *****************************************************************************

 Declaration of interfaces for:
 - a R1CS constraint,
 - a R1CS variable assignment, and
 - a R1CS constraint system.

 Above, R1CS stands for "Rank-1 Constraint System".

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef R1CS_HPP_
#define R1CS_HPP_

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <ringsnark/relations/variable.hpp>

namespace ringsnark {

/************************* R1CS constraint ***********************************/

    template<typename RingT>
    class r1cs_constraint;

    template<typename RingT>
    std::ostream &operator<<(std::ostream &out, const r1cs_constraint<RingT> &c);

    template<typename RingT>
    std::istream &operator>>(std::istream &in, r1cs_constraint<RingT> &c);

/**
 * A R1CS constraint is a formal expression of the form
 *
 *                < A , X > * < B , X > = < C , X > ,
 *
 * where X = (x_0,x_1,...,x_m) is a vector of formal variables and A,B,C each
 * consist of 1+m elements in <RingT>.
 *
 * A R1CS constraint is used to construct a R1CS constraint system (see below).
 */
    template<typename RingT>
    class r1cs_constraint {
    public:

        linear_combination<RingT> a, b, c;

        r1cs_constraint(const r1cs_constraint<RingT> &other) = default;

        r1cs_constraint() = default;

        r1cs_constraint(const linear_combination<RingT> &a,
                        const linear_combination<RingT> &b,
                        const linear_combination<RingT> &c);

        r1cs_constraint(const std::initializer_list<linear_combination<RingT> > &A,
                        const std::initializer_list<linear_combination<RingT> > &B,
                        const std::initializer_list<linear_combination<RingT> > &C);

        bool operator==(const r1cs_constraint<RingT> &other) const;

        friend std::ostream &operator<<<RingT>(std::ostream &out, const r1cs_constraint<RingT> &c);

        friend std::istream &operator>><RingT>(std::istream &in, r1cs_constraint<RingT> &c);
    };

/************************* R1CS variable assignment **************************/

/**
 * A R1CS variable assignment is a vector of <RingT> elements that represents
 * a candidate solution to a R1CS constraint system (see below).
 * This does not include the constant input 1, which is handled separately.
 */

    template<typename RingT>
    using r1cs_primary_input = std::vector<RingT>;

    template<typename RingT>
    using r1cs_auxiliary_input = std::vector<RingT>;

    template<typename RingT>
    using r1cs_variable_assignment = std::vector<RingT>;

/************************* R1CS constraint system ****************************/

    template<typename RingT>
    class r1cs_constraint_system;

    template<typename RingT>
    std::ostream &operator<<(std::ostream &out, const r1cs_constraint_system<RingT> &cs);

    template<typename RingT>
    std::istream &operator>>(std::istream &in, r1cs_constraint_system<RingT> &cs);

/**
 * A system of R1CS constraints looks like
 *
 *     { < A_k , X > * < B_k , X > = < C_k , X > }_{k=1}^{n}  .
 *
 * In other words, the system is satisfied if and only if there exist a
 * USCS variable assignment for which each R1CS constraint is satisfied.
 *
 * NOTE:
 * The 0-th variable (i.e., "x_{0}") always represents the constant 1.
 * Thus, the 0-th variable is not included in num_variables.
 */
    template<typename RingT>
    class r1cs_constraint_system {
    public:
        size_t primary_input_size;
        size_t auxiliary_input_size;

        std::vector<r1cs_constraint<RingT> > constraints;

        r1cs_constraint_system(const r1cs_constraint_system<RingT> &other) = default;

        r1cs_constraint_system() = default;

        [[nodiscard]] size_t num_inputs() const;

        [[nodiscard]] size_t num_variables() const;

        [[nodiscard]] size_t num_constraints() const;

#ifdef DEBUG
        std::map<size_t, std::string> constraint_annotations;
        std::map<size_t, std::string> variable_annotations;
#endif

        [[nodiscard]] bool is_valid() const;

        [[nodiscard]] bool is_satisfied(const r1cs_primary_input<RingT> &primary_input,
                                        const r1cs_auxiliary_input<RingT> &auxiliary_input) const;

        void add_constraint(const r1cs_constraint<RingT> &c);

        void add_constraint(const r1cs_constraint<RingT> &c, const std::string &annotation);

        void swap_AB_if_beneficial();

        bool operator==(const r1cs_constraint_system<RingT> &other) const;

        friend std::ostream &operator
        <<<RingT>(
        std::ostream &out,
        const r1cs_constraint_system<RingT> &cs
        );

        friend std::istream &operator>><RingT>(std::istream &in, r1cs_constraint_system<RingT> &cs);

        void report_linear_constraint_statistics() const;
    };


} // ringsnark

#include <ringsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.tcc>

#endif // R1CS_HPP_
