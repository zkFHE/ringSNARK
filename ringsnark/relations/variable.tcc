/** @file
 *****************************************************************************
 Implementation of interfaces for:
 - a variable (i.e., x_i),
 - a linear term (i.e., a_i * x_i), and
 - a linear combination (i.e., sum_i a_i * x_i).
 See variable.hpp .
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef VARIABLE_TCC_
#define VARIABLE_TCC_

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "variable.hpp"

namespace ringsnark {

    template<typename RingT>
    linear_term<RingT> variable<RingT>::operator*(const integer_coeff_t int_coeff) const {
        return linear_term<RingT>(*this, int_coeff);
    }

    template<typename RingT>
    linear_term<RingT> variable<RingT>::operator*(const RingT &field_coeff) const {
        return linear_term<RingT>(*this, field_coeff);
    }

    template<typename RingT>
    linear_combination<RingT> variable<RingT>::operator+(const linear_combination<RingT> &other) const {
        linear_combination<RingT> result;

        result.add_term(*this);
        result.terms.insert(result.terms.begin(), other.terms.begin(), other.terms.end());

        return result;
    }

    template<typename RingT>
    linear_combination<RingT> variable<RingT>::operator-(const linear_combination<RingT> &other) const {
        return (*this) + (-other);
    }

    template<typename RingT>
    linear_term<RingT> variable<RingT>::operator-() const {
        return linear_term<RingT>(*this, -RingT::one());
    }

    template<typename RingT>
    bool variable<RingT>::operator==(const variable<RingT> &other) const {
        return (this->index == other.index);
    }

    template<typename RingT>
    linear_term<RingT> operator*(const integer_coeff_t int_coeff, const variable<RingT> &var) {
        return linear_term<RingT>(var, int_coeff);
    }

    template<typename RingT>
    linear_term<RingT> operator*(const RingT &field_coeff, const variable<RingT> &var) {
        return linear_term<RingT>(var, field_coeff);
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const variable<RingT> &var) {
        return linear_combination<RingT>(int_coeff) + var;
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &field_coeff, const variable<RingT> &var) {
        return linear_combination<RingT>(field_coeff) + var;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const variable<RingT> &var) {
        return linear_combination<RingT>(int_coeff) - var;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &field_coeff, const variable<RingT> &var) {
        return linear_combination<RingT>(field_coeff) - var;
    }

    template<typename RingT>
    linear_term<RingT>::linear_term(const variable<RingT> &var) :
            index(var.index), coeff(RingT::one()) {
    }

    template<typename RingT>
    linear_term<RingT>::linear_term(const variable<RingT> &var, const integer_coeff_t int_coeff) :
            index(var.index), coeff(RingT(int_coeff)) {
    }

    template<typename RingT>
    linear_term<RingT>::linear_term(const variable<RingT> &var, const RingT &coeff) :
            index(var.index), coeff(coeff) {
    }

    template<typename RingT>
    linear_term<RingT> linear_term<RingT>::operator*(const integer_coeff_t int_coeff) const {
        return (this->operator*(RingT(int_coeff)));
    }

    template<typename RingT>
    linear_term<RingT> linear_term<RingT>::operator*(const RingT &field_coeff) const {
        return linear_term<RingT>(this->index, field_coeff * this->coeff);
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const linear_term<RingT> &lt) {
        return linear_combination<RingT>(int_coeff) + lt;
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &field_coeff, const linear_term<RingT> &lt) {
        return linear_combination<RingT>(field_coeff) + lt;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const linear_term<RingT> &lt) {
        return linear_combination<RingT>(int_coeff) - lt;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &field_coeff, const linear_term<RingT> &lt) {
        return linear_combination<RingT>(field_coeff) - lt;
    }

    template<typename RingT>
    linear_combination<RingT> linear_term<RingT>::operator+(const linear_combination<RingT> &other) const {
        return linear_combination<RingT>(*this) + other;
    }

    template<typename RingT>
    linear_combination<RingT> linear_term<RingT>::operator-(const linear_combination<RingT> &other) const {
        return (*this) + (-other);
    }

    template<typename RingT>
    linear_term<RingT> linear_term<RingT>::operator-() const {
        return linear_term<RingT>(this->index, -this->coeff);
    }

    template<typename RingT>
    bool linear_term<RingT>::operator==(const linear_term<RingT> &other) const {
        return (this->index == other.index &&
                this->coeff == other.coeff);
    }

    template<typename RingT>
    linear_term<RingT> operator*(const integer_coeff_t int_coeff, const linear_term<RingT> &lt) {
        return RingT(int_coeff) * lt;
    }

    template<typename RingT>
    linear_term<RingT> operator*(const RingT &field_coeff, const linear_term<RingT> &lt) {
        return linear_term<RingT>(lt.index, field_coeff * lt.coeff);
    }

    template<typename RingT>
    linear_combination<RingT>::linear_combination(const integer_coeff_t int_coeff) {
        this->add_term(linear_term<RingT>(0, int_coeff));
    }

    template<typename RingT>
    linear_combination<RingT>::linear_combination(const RingT &field_coeff) {
        this->add_term(linear_term<RingT>(0, field_coeff));
    }

    template<typename RingT>
    linear_combination<RingT>::linear_combination(const variable<RingT> &var) {
        this->add_term(var);
    }

    template<typename RingT>
    linear_combination<RingT>::linear_combination(const linear_term<RingT> &lt) {
        this->add_term(lt);
    }

    template<typename RingT>
    typename std::vector<linear_term<RingT> >

    ::const_iterator linear_combination<RingT>::begin() const {
        return terms.begin();
    }

    template<typename RingT>
    typename std::vector<linear_term<RingT> >

    ::const_iterator linear_combination<RingT>::end() const {
        return terms.end();
    }

    template<typename RingT>
    void linear_combination<RingT>::add_term(const variable<RingT> &var) {
        this->terms.emplace_back(linear_term<RingT>(var.index, RingT::one()));
    }

    template<typename RingT>
    void linear_combination<RingT>::add_term(const variable<RingT> &var, const integer_coeff_t int_coeff) {
        this->terms.emplace_back(linear_term<RingT>(var.index, int_coeff));
    }

    template<typename RingT>
    void linear_combination<RingT>::add_term(const variable<RingT> &var, const RingT &coeff) {
        this->terms.emplace_back(linear_term<RingT>(var.index, coeff));
    }

    template<typename RingT>
    void linear_combination<RingT>::add_term(const linear_term<RingT> &other) {
        this->terms.emplace_back(other);
    }

    template<typename RingT>
    linear_combination<RingT> linear_combination<RingT>::operator*(const integer_coeff_t int_coeff) const {
        return (*this) * RingT(int_coeff);
    }

    template<typename RingT>
    RingT linear_combination<RingT>::evaluate(const std::vector<RingT> &assignment) const {
        RingT acc = RingT::zero();
        for (auto &lt: terms) {
            acc += (lt.index == 0 ? RingT::one() : assignment.at(lt.index - 1)) * lt.coeff;
        }
        return acc;
    }

    template<typename RingT>
    linear_combination<RingT> linear_combination<RingT>::operator*(const RingT &field_coeff) const {
        linear_combination<RingT> result;
        result.terms.reserve(this->terms.size());
        for (const linear_term<RingT> &lt: this->terms) {
            result.terms.emplace_back(lt * field_coeff);
        }
        return result;
    }

    template<typename RingT>
    linear_combination<RingT> linear_combination<RingT>::operator+(const linear_combination<RingT> &other) const {
        linear_combination<RingT> result;

        auto it1 = this->terms.begin();
        auto it2 = other.terms.begin();

        /* invariant: it1 and it2 always point to unprocessed items in the corresponding linear combinations */
        while (it1 != this->terms.end() && it2 != other.terms.end()) {
            if (it1->index < it2->index) {
                result.terms.emplace_back(*it1);
                ++it1;
            } else if (it1->index > it2->index) {
                result.terms.emplace_back(*it2);
                ++it2;
            } else {
                /* it1->index == it2->index */
                result.terms.emplace_back(linear_term<RingT>(variable<RingT>(it1->index), it1->coeff + it2->coeff));
                ++it1;
                ++it2;
            }
        }

        if (it1 != this->terms.end()) {
            result.terms.insert(result.terms.end(), it1, this->terms.end());
        } else {
            result.terms.insert(result.terms.end(), it2, other.terms.end());
        }

        return result;
    }

    template<typename RingT>
    linear_combination<RingT> linear_combination<RingT>::operator-(const linear_combination<RingT> &other) const {
        return (*this) + (-other);
    }

    template<typename RingT>
    linear_combination<RingT> linear_combination<RingT>::operator-() const {
        return (*this) * (-RingT::one());
    }

    template<typename RingT>
    bool linear_combination<RingT>::operator==(const linear_combination<RingT> &other) const {
        return (this->terms == other.terms);
    }

    template<typename RingT>
    bool linear_combination<RingT>::is_valid(const size_t num_variables) const {
        /* check that all terms in linear combination are sorted */
        for (size_t i = 1; i < terms.size(); ++i) {
            if (terms[i - 1].index >= terms[i].index) {
                return false;
            }
        }

        /* check that the variables are in proper range. as the variables
           are sorted, it suffices to check the last term */
        if ((--terms.end())->index >= num_variables) {
            return false;
        }

        return true;
    }

    template<typename RingT>
    void linear_combination<RingT>::print(const std::map<size_t, std::string> &variable_annotations) const {
        for (auto &lt: terms) {
            if (lt.index == 0) {
                printf("    1 * ");
                lt.coeff.print();
            } else {
                auto it = variable_annotations.find(lt.index);
                printf("    x_%zu (%s) * ", lt.index,
                       (it == variable_annotations.end() ? "no annotation" : it->second.c_str()));
                lt.coeff.print();
            }
        }
    }

    template<typename RingT>
    void linear_combination<RingT>::print_with_assignment(const std::vector<RingT> &full_assignment,
                                                          const std::map<size_t, std::string> &variable_annotations) const {
        for (auto &lt: terms) {
            if (lt.index == 0) {
                printf("    1 * ");
                lt.coeff.print();
            } else {
                printf("    x_%zu * ", lt.index);
                lt.coeff.print();

                auto it = variable_annotations.find(lt.index);
                printf("    where x_%zu (%s) was assigned value ", lt.index,
                       (it == variable_annotations.end() ? "no annotation" : it->second.c_str()));
                full_assignment[lt.index - 1].print();
                printf("      i.e. negative of ");
                (-full_assignment[lt.index - 1]).print();
            }
        }
    }

    template<typename RingT>
    std::ostream &operator<<(std::ostream &out, const linear_combination<RingT> &lc) {
        auto OUTPUT_NEWLINE = "\n";
        out << lc.terms.size() << "\n";
        for (const linear_term<RingT> &lt: lc.terms) {
            out << lt.index << "\n";
            out << lt.coeff << OUTPUT_NEWLINE;
        }

        return out;
    }

    template<typename RingT>
    std::istream &operator>>(std::istream &in, linear_combination<RingT> &lc) {
        /*
        lc.terms.clear();

        size_t s;
        in >> s;

        libff::consume_newline(in);

        lc.terms.reserve(s);

        for (size_t i = 0; i < s; ++i) {
            linear_term<RingT> lt;
            in >> lt.index;
            libff::consume_newline(in);
            in >> lt.coeff;
            libff::consume_OUTPUT_NEWLINE(in);
            lc.terms.emplace_back(lt);
        }

        return in;
         */
        throw std::runtime_error("Not Implemented");
    }

    template<typename RingT>
    linear_combination<RingT> operator*(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc) {
        return lc * int_coeff;
    }

    template<typename RingT>
    linear_combination<RingT> operator*(const RingT &field_coeff, const linear_combination<RingT> &lc) {
        return lc * field_coeff;
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc) {
        return linear_combination<RingT>(int_coeff) + lc;
    }

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &field_coeff, const linear_combination<RingT> &lc) {
        return linear_combination<RingT>(field_coeff) + lc;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc) {
        return linear_combination<RingT>(int_coeff) - lc;
    }

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &field_coeff, const linear_combination<RingT> &lc) {
        return linear_combination<RingT>(field_coeff) - lc;
    }

    template<typename RingT>
    linear_combination<RingT>::linear_combination(const std::vector<linear_term<RingT>> &all_terms) {
        if (all_terms.empty()) {
            return;
        }

        terms = all_terms;
        std::sort(terms.begin(), terms.end(),
                  [](linear_term<RingT> a, linear_term<RingT> b) { return a.index < b.index; });

        auto result_it = terms.begin();
        for (auto it = ++terms.begin(); it != terms.end(); ++it) {
            if (it->index == result_it->index) {
                result_it->coeff += it->coeff;
            } else {
                *(++result_it) = *it;
            }
        }
        terms.resize((result_it - terms.begin()) + 1);
    }

} // ringsnark

#endif // VARIABLE_TCC

