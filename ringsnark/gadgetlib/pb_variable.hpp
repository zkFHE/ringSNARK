/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PB_VARIABLE_HPP_
#define PB_VARIABLE_HPP_

#include <cstddef>
#include <string>
#include <vector>

#include <ringsnark/relations/variable.hpp>

namespace ringsnark {

typedef size_t lc_index_t;

template<typename RingT>
class protoboard;

template<typename RingT>
class pb_variable : public variable<RingT> {
public:
    pb_variable(const var_index_t index = 0) : variable<RingT>(index) {};

    void allocate(protoboard<RingT> &pb, const std::string &annotation="");
};

template<typename RingT>
class pb_variable_array : private std::vector<pb_variable<RingT> >
{
    typedef std::vector<pb_variable<RingT> > contents;
public:
    using typename contents::iterator;
    using typename contents::const_iterator;
    using typename contents::reverse_iterator;
    using typename contents::const_reverse_iterator;

    using contents::begin;
    using contents::end;
    using contents::rbegin;
    using contents::rend;
    using contents::emplace_back;
    using contents::insert;
    using contents::reserve;
    using contents::size;
    using contents::empty;
    using contents::operator[];
    using contents::resize;

    pb_variable_array() : contents() {};
    pb_variable_array(size_t count, const pb_variable<RingT> &value) : contents(count, value) {};
    pb_variable_array(typename contents::const_iterator first, typename contents::const_iterator last) : contents(first, last) {};
    pb_variable_array(typename contents::const_reverse_iterator first, typename contents::const_reverse_iterator last) : contents(first, last) {};
    void allocate(protoboard<RingT> &pb, const size_t n, const std::string &annotation_prefix="");

    void fill_with_field_elements(protoboard<RingT> &pb, const std::vector<RingT>& vals) const;
//    void fill_with_bits(protoboard<RingT> &pb, const libff::bit_vector& bits) const;
    void fill_with_bits_of_ulong(protoboard<RingT> &pb, const unsigned long i) const;
    void fill_with_bits_of_field_element(protoboard<RingT> &pb, const RingT &r) const;

    std::vector<RingT> get_vals(const protoboard<RingT> &pb) const;
//    libff::bit_vector get_bits(const protoboard<RingT> &pb) const;

    RingT get_field_element_from_bits(const protoboard<RingT> &pb) const;
};

/* index 0 corresponds to the constant term (used in legacy code) */
#define ONE pb_variable<RingT>(0)

template<typename RingT>
class pb_linear_combination : public linear_combination<RingT> {
public:
    bool is_variable;
    lc_index_t index;

    pb_linear_combination();
    pb_linear_combination(const pb_variable<RingT> &var);

    void assign(protoboard<RingT> &pb, const linear_combination<RingT> &lc);
    void evaluate(protoboard<RingT> &pb) const;

    bool is_constant() const;
    RingT constant_term() const;
};

template<typename RingT>
class pb_linear_combination_array : private std::vector<pb_linear_combination<RingT> >
{
    typedef std::vector<pb_linear_combination<RingT> > contents;
public:
    using typename contents::iterator;
    using typename contents::const_iterator;
    using typename contents::reverse_iterator;
    using typename contents::const_reverse_iterator;

    using contents::begin;
    using contents::end;
    using contents::rbegin;
    using contents::rend;
    using contents::emplace_back;
    using contents::insert;
    using contents::reserve;
    using contents::size;
    using contents::empty;
    using contents::operator[];
    using contents::resize;

    pb_linear_combination_array() : contents() {};
    pb_linear_combination_array(const pb_variable_array<RingT> &arr) { for (auto &v : arr) this->emplace_back(pb_linear_combination<RingT>(v)); };
    pb_linear_combination_array(size_t count) : contents(count) {};
    pb_linear_combination_array(size_t count, const pb_linear_combination<RingT> &value) : contents(count, value) {};
    pb_linear_combination_array(typename contents::const_iterator first, typename contents::const_iterator last) : contents(first, last) {};
    pb_linear_combination_array(typename contents::const_reverse_iterator first, typename contents::const_reverse_iterator last) : contents(first, last) {};

    void evaluate(protoboard<RingT> &pb) const;

    void fill_with_field_elements(protoboard<RingT> &pb, const std::vector<RingT>& vals) const;
//    void fill_with_bits(protoboard<RingT> &pb, const libff::bit_vector& bits) const;
    void fill_with_bits_of_ulong(protoboard<RingT> &pb, const unsigned long i) const;
    void fill_with_bits_of_field_element(protoboard<RingT> &pb, const RingT &r) const;

    std::vector<RingT> get_vals(const protoboard<RingT> &pb) const;
//    libff::bit_vector get_bits(const protoboard<RingT> &pb) const;

    RingT get_field_element_from_bits(const protoboard<RingT> &pb) const;
};

template<typename RingT>
linear_combination<RingT> pb_sum(const pb_linear_combination_array<RingT> &v);

template<typename RingT>
linear_combination<RingT> pb_packing_sum(const pb_linear_combination_array<RingT> &v);

template<typename RingT>
linear_combination<RingT> pb_coeff_sum(const pb_linear_combination_array<RingT> &v, const std::vector<RingT> &coeffs);

} // ringsnark
#include <ringsnark/gadgetlib/pb_variable.tcc>

#endif // PB_VARIABLE_HPP_
