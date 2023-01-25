/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PROTOBOARD_HPP_
#define PROTOBOARD_HPP_

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <ringsnark/gadgetlib/pb_variable.hpp>
#include <ringsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>

namespace ringsnark {

    template<typename RingT>
    class r1cs_constraint;

    template<typename RingT>
    class r1cs_constraint_system;

    template<typename RingT>
    class protoboard {
    private:
        RingT constant_term; /* only here, because pb.val() needs to be able to return reference to the constant 1 term */
        r1cs_variable_assignment<RingT> values; /* values[0] will hold the value of the first allocated variable of the protoboard, *NOT* constant 1 */
        var_index_t next_free_var;
        lc_index_t next_free_lc;
        std::vector<RingT> lc_values;
        r1cs_constraint_system<RingT> constraint_system;

    public:
        protoboard();

        void clear_values();

        RingT &val(const pb_variable<RingT> &var);

        RingT val(const pb_variable<RingT> &var) const;

        RingT &lc_val(const pb_linear_combination<RingT> &lc);

        RingT lc_val(const pb_linear_combination<RingT> &lc) const;

        void add_r1cs_constraint(const r1cs_constraint<RingT> &constr, const std::string &annotation = "");

        void augment_variable_annotation(const pb_variable<RingT> &v, const std::string &postfix);

        bool is_satisfied() const;

        void dump_variables() const;

        size_t num_constraints() const;

        size_t num_inputs() const;

        size_t num_variables() const;

        void set_input_sizes(const size_t primary_input_size);

        r1cs_variable_assignment<RingT> full_variable_assignment() const;

        r1cs_primary_input<RingT> primary_input() const;

        r1cs_auxiliary_input<RingT> auxiliary_input() const;

        r1cs_constraint_system<RingT> get_constraint_system() const;

        friend class pb_variable<RingT>;

        friend class pb_linear_combination<RingT>;

    private:
        var_index_t allocate_var_index(const std::string &annotation = "");

        lc_index_t allocate_lc_index();
    };

} // ringsnark
#include <ringsnark/gadgetlib/protoboard.tcc>

#endif // PROTOBOARD_HPP_
