/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PB_VARIABLE_TCC_
#define PB_VARIABLE_TCC_

#include <cassert>

#include <ringsnark/gadgetlib/protoboard.hpp>

namespace ringsnark {

    template<typename RingT>
    void pb_variable<RingT>::allocate(protoboard<RingT> &pb, const std::string &annotation) {
        this->index = pb.allocate_var_index(annotation);
    }

/* allocates pb_variable<RingT> array in MSB->LSB order */
    template<typename RingT>
    void
    pb_variable_array<RingT>::allocate(protoboard<RingT> &pb, const size_t n, const std::string &annotation_prefix) {
#ifdef DEBUG
        assert(annotation_prefix != "");
#endif
        (*this).resize(n);

        for (size_t i = 0; i < n; ++i) {
//        (*this)[i].allocate(pb, FMT(annotation_prefix, "_%zu", i));
            (*this)[i].allocate(pb, annotation_prefix + std::to_string(i));
        }
    }

    template<typename RingT>
    void
    pb_variable_array<RingT>::fill_with_field_elements(protoboard<RingT> &pb, const std::vector<RingT> &vals) const {
        assert(this->size() == vals.size());
        for (size_t i = 0; i < vals.size(); ++i) {
            pb.val((*this)[i]) = vals[i];
        }
    }

//template<typename RingT>
//void pb_variable_array<RingT>::fill_with_bits(protoboard<RingT> &pb, const libff::bit_vector& bits) const
//{
//    assert(this->size() == bits.size());
//    for (size_t i = 0; i < bits.size(); ++i)
//    {
//        pb.val((*this)[i]) = (bits[i] ? RingT::one() : RingT::zero());
//    }
//}

    template<typename RingT>
    void pb_variable_array<RingT>::fill_with_bits_of_field_element(protoboard<RingT> &pb, const RingT &r) const {
//    const libff::bigint<RingT::num_limbs> rint = r.as_bigint();
//    for (size_t i = 0; i < this->size(); ++i)
//    {
//        pb.val((*this)[i]) = rint.test_bit(i) ? RingT::one() : RingT::zero();
//    }
    }

    template<typename RingT>
    void pb_variable_array<RingT>::fill_with_bits_of_ulong(protoboard<RingT> &pb, const unsigned long i) const {
        this->fill_with_bits_of_field_element(pb, RingT(i, true));
    }

    template<typename RingT>
    std::vector<RingT> pb_variable_array<RingT>::get_vals(const protoboard<RingT> &pb) const {
        std::vector<RingT> result(this->size());
        for (size_t i = 0; i < this->size(); ++i) {
            result[i] = pb.val((*this)[i]);
        }
        return result;
    }
//
//template<typename RingT>
//libff::bit_vector pb_variable_array<RingT>::get_bits(const protoboard<RingT> &pb) const
//{
//    libff::bit_vector result;
//    for (size_t i = 0; i < this->size(); ++i)
//    {
//        const RingT v = pb.val((*this)[i]);
//        assert(v == RingT::zero() || v == RingT::one());
//        result.push_back(v == RingT::one());
//    }
//    return result;
//}

    template<typename RingT>
    RingT pb_variable_array<RingT>::get_field_element_from_bits(const protoboard<RingT> &pb) const {
        RingT result = RingT::zero();

        for (size_t i = 0; i < this->size(); ++i) {
            /* push in the new bit */
            const RingT v = pb.val((*this)[this->size() - 1 - i]);
            assert(v == RingT::zero() || v == RingT::one());
            result += result + v;
        }

        return result;
    }

    template<typename RingT>
    pb_linear_combination<RingT>::pb_linear_combination() {
        this->is_variable = false;
    }

    template<typename RingT>
    pb_linear_combination<RingT>::pb_linear_combination(const pb_variable<RingT> &var) {
        this->is_variable = true;
        this->index = var.index;
        this->terms.emplace_back(linear_term<RingT>(var));
    }

    template<typename RingT>
    void pb_linear_combination<RingT>::assign(protoboard<RingT> &pb, const linear_combination<RingT> &lc) {
        assert(this->is_variable == false);
        this->index = pb.allocate_lc_index();
        this->terms = lc.terms;
    }

    template<typename RingT>
    void pb_linear_combination<RingT>::evaluate(protoboard<RingT> &pb) const {
        if (this->is_variable) {
            return; // do nothing
        }

        RingT sum = 0;
        for (auto term: this->terms) {
            sum += term.coeff * pb.val(pb_variable<RingT>(term.index));
        }

        pb.lc_val(*this) = sum;
    }

    template<typename RingT>
    bool pb_linear_combination<RingT>::is_constant() const {
        if (is_variable) {
            return (index == 0);
        } else {
            for (auto term: this->terms) {
                if (term.index != 0) {
                    return false;
                }
            }

            return true;
        }
    }

    template<typename RingT>
    RingT pb_linear_combination<RingT>::constant_term() const {
        if (is_variable) {
            return (index == 0 ? RingT::one() : RingT::zero());
        } else {
            RingT result = RingT::zero();
            for (auto term: this->terms) {
                if (term.index == 0) {
                    result += term.coeff;
                }
            }
            return result;
        }
    }

    template<typename RingT>
    void pb_linear_combination_array<RingT>::evaluate(protoboard<RingT> &pb) const {
        for (size_t i = 0; i < this->size(); ++i) {
            (*this)[i].evaluate(pb);
        }
    }

    template<typename RingT>
    void pb_linear_combination_array<RingT>::fill_with_field_elements(protoboard<RingT> &pb,
                                                                      const std::vector<RingT> &vals) const {
        assert(this->size() == vals.size());
        for (size_t i = 0; i < vals.size(); ++i) {
            pb.lc_val((*this)[i]) = vals[i];
        }
    }

//template<typename RingT>
//void pb_linear_combination_array<RingT>::fill_with_bits(protoboard<RingT> &pb, const libff::bit_vector& bits) const
//{
//    assert(this->size() == bits.size());
//    for (size_t i = 0; i < bits.size(); ++i)
//    {
//        pb.lc_val((*this)[i]) = (bits[i] ? RingT::one() : RingT::zero());
//    }
//}

    template<typename RingT>
    void
    pb_linear_combination_array<RingT>::fill_with_bits_of_field_element(protoboard<RingT> &pb, const RingT &r) const {
//    const libff::bigint<RingT::num_limbs> rint = r.as_bigint();
//    for (size_t i = 0; i < this->size(); ++i)
//    {
//        pb.lc_val((*this)[i]) = rint.test_bit(i) ? RingT::one() : RingT::zero();
//    }
    }

    template<typename RingT>
    void
    pb_linear_combination_array<RingT>::fill_with_bits_of_ulong(protoboard<RingT> &pb, const unsigned long i) const {
        this->fill_with_bits_of_field_element(pb, RingT(i));
    }

    template<typename RingT>
    std::vector<RingT> pb_linear_combination_array<RingT>::get_vals(const protoboard<RingT> &pb) const {
        std::vector<RingT> result(this->size());
        for (size_t i = 0; i < this->size(); ++i) {
            result[i] = pb.lc_val((*this)[i]);
        }
        return result;
    }

//template<typename RingT>
//libff::bit_vector pb_linear_combination_array<RingT>::get_bits(const protoboard<RingT> &pb) const
//{
//    libff::bit_vector result;
//    for (size_t i = 0; i < this->size(); ++i)
//    {
//        const RingT v = pb.lc_val((*this)[i]);
//        assert(v == RingT::zero() || v == RingT::one());
//        result.push_back(v == RingT::one());
//    }
//    return result;
//}

    template<typename RingT>
    RingT pb_linear_combination_array<RingT>::get_field_element_from_bits(const protoboard<RingT> &pb) const {
        RingT result = RingT::zero();

        for (size_t i = 0; i < this->size(); ++i) {
            /* push in the new bit */
            const RingT v = pb.lc_val((*this)[this->size() - 1 - i]);
            assert(v == RingT::zero() || v == RingT::one());
            result += result + v;
        }

        return result;
    }

    template<typename RingT>
    linear_combination<RingT> pb_sum(const pb_linear_combination_array<RingT> &v) {
        linear_combination<RingT> result;
        for (auto &term: v) {
            result = result + term;
        }

        return result;
    }

    template<typename RingT>
    linear_combination<RingT> pb_packing_sum(const pb_linear_combination_array<RingT> &v) {
        RingT twoi = RingT::one(); // will hold 2^i entering each iteration
        std::vector<linear_term<RingT> > all_terms;
        for (auto &lc: v) {
            for (auto &term: lc.terms) {
                all_terms.emplace_back(twoi * term);
            }
            twoi += twoi;
        }

        return linear_combination<RingT>(all_terms);
    }

    template<typename RingT>
    linear_combination<RingT>
    pb_coeff_sum(const pb_linear_combination_array<RingT> &v, const std::vector<RingT> &coeffs) {
        assert(v.size() == coeffs.size());
        std::vector<linear_term<RingT> > all_terms;

        auto coeff_it = coeffs.begin();
        for (auto &lc: v) {
            for (auto &term: lc.terms) {
                all_terms.emplace_back((*coeff_it) * term);
            }
            ++coeff_it;
        }

        return linear_combination<RingT>(all_terms);
    }


} // ringsnark
#endif // PB_VARIABLE_TCC
