/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PROTOBOARD_TCC_
#define PROTOBOARD_TCC_

#include <cstdarg>
#include <cstdio>

namespace ringsnark {

    template<typename RingT>
    protoboard<RingT>::protoboard() {
        constant_term = RingT::one();

#ifdef DEBUG
        constraint_system.variable_annotations[0] = "ONE";
#endif

        next_free_var = 1; /* to account for constant 1 term */
        next_free_lc = 0;
    }

    template<typename RingT>
    void protoboard<RingT>::clear_values() {
        std::fill(values.begin(), values.end(), RingT::zero());
    }

    template<typename RingT>
    var_index_t protoboard<RingT>::allocate_var_index(const std::string &annotation) {
#ifdef DEBUG
        assert(annotation != "");
        constraint_system.variable_annotations[next_free_var] = annotation;
#else
//    libff::UNUSED(annotation);
#endif
        ++constraint_system.auxiliary_input_size;
        values.emplace_back(RingT::zero());
        return next_free_var++;
    }

    template<typename RingT>
    lc_index_t protoboard<RingT>::allocate_lc_index() {
        lc_values.emplace_back(RingT::zero());
        return next_free_lc++;
    }

    template<typename RingT>
    RingT &protoboard<RingT>::val(const pb_variable <RingT> &var) {
        assert(var.index <= values.size());
        return (var.index == 0 ? constant_term : values[var.index - 1]);
    }

    template<typename RingT>
    RingT protoboard<RingT>::val(const pb_variable <RingT> &var) const {
        assert(var.index <= values.size());
        return (var.index == 0 ? constant_term : values[var.index - 1]);
    }

    template<typename RingT>
    RingT &protoboard<RingT>::lc_val(const pb_linear_combination <RingT> &lc) {
        if (lc.is_variable) {
            return this->val(pb_variable<RingT>(lc.index));
        } else {
            assert(lc.index < lc_values.size());
            return lc_values[lc.index];
        }
    }

    template<typename RingT>
    RingT protoboard<RingT>::lc_val(const pb_linear_combination <RingT> &lc) const {
        if (lc.is_variable) {
            return this->val(pb_variable<RingT>(lc.index));
        } else {
            assert(lc.index < lc_values.size());
            return lc_values[lc.index];
        }
    }

    template<typename RingT>
    void protoboard<RingT>::add_r1cs_constraint(const r1cs_constraint <RingT> &constr, const std::string &annotation) {
#ifdef DEBUG
        assert(annotation != "");
        constraint_system.constraint_annotations[constraint_system.constraints.size()] = annotation;
#else
//    libff::UNUSED(annotation);
#endif
        constraint_system.constraints.emplace_back(constr);
    }

    template<typename RingT>
    void protoboard<RingT>::augment_variable_annotation(const pb_variable <RingT> &v, const std::string &postfix) {
#ifdef DEBUG
        auto it = constraint_system.variable_annotations.find(v.index);
        constraint_system.variable_annotations[v.index] = (it == constraint_system.variable_annotations.end() ? "" : it->second + " ") + postfix;
#endif
    }

    template<typename RingT>
    bool protoboard<RingT>::is_satisfied() const {
        return constraint_system.is_satisfied(primary_input(), auxiliary_input());
    }

    template<typename RingT>
    void protoboard<RingT>::dump_variables() const {
#ifdef DEBUG
        for (size_t i = 0; i < constraint_system.num_variables; ++i)
        {
            printf("%-40s --> ", constraint_system.variable_annotations[i].c_str());
            values[i].as_bigint().print_hex();
        }
#endif
    }

    template<typename RingT>
    size_t protoboard<RingT>::num_constraints() const {
        return constraint_system.num_constraints();
    }

    template<typename RingT>
    size_t protoboard<RingT>::num_inputs() const {
        return constraint_system.num_inputs();
    }

    template<typename RingT>
    size_t protoboard<RingT>::num_variables() const {
        return next_free_var - 1;
    }

    template<typename RingT>
    void protoboard<RingT>::set_input_sizes(const size_t primary_input_size) {
        assert(primary_input_size <= num_variables());
        constraint_system.primary_input_size = primary_input_size;
        constraint_system.auxiliary_input_size = num_variables() - primary_input_size;
    }

    template<typename RingT>
    r1cs_variable_assignment <RingT> protoboard<RingT>::full_variable_assignment() const {
        return values;
    }

    template<typename RingT>
    r1cs_primary_input <RingT> protoboard<RingT>::primary_input() const {
        return r1cs_primary_input<RingT>(values.begin(), values.begin() + num_inputs());
    }

    template<typename RingT>
    r1cs_auxiliary_input <RingT> protoboard<RingT>::auxiliary_input() const {
        return r1cs_auxiliary_input<RingT>(values.begin() + num_inputs(), values.end());
    }

    template<typename RingT>
    r1cs_constraint_system <RingT> protoboard<RingT>::get_constraint_system() const {
        return constraint_system;
    }

} // ringsnark
#endif // PROTOBOARD_TCC_
