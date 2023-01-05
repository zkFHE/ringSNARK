/** @file
 *****************************************************************************

 Declaration of interfaces for:
 - a variable (i.e., x_i),
 - a linear term (i.e., a_i * x_i), and
 - a linear combination (i.e., sum_i a_i * x_i).

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef VARIABLE_HPP_
#define VARIABLE_HPP_

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace ringsnark {

/**
 * Mnemonic typedefs.
 */
    typedef size_t var_index_t;
    typedef long integer_coeff_t;

/**
 * Forward declaration.
 */
    template<typename RingT>
    class linear_term;

/**
 * Forward declaration.
 */
    template<typename RingT>
    class linear_combination;

/********************************* Variable **********************************/

/**
 * A variable represents a formal expression of the form "x_{index}".
 */
    template<typename RingT>
    class variable {
    public:

        var_index_t index;

        variable(const var_index_t index = 0) : index(index) {};

        linear_term<RingT> operator*(const integer_coeff_t int_coeff) const;

        linear_term<RingT> operator*(const RingT &ring_coeff) const;

        linear_combination<RingT> operator+(const linear_combination<RingT> &other) const;

        linear_combination<RingT> operator-(const linear_combination<RingT> &other) const;

        linear_term<RingT> operator-() const;

        bool operator==(const variable<RingT> &other) const;
    };

    template<typename RingT>
    linear_term<RingT> operator*(const integer_coeff_t int_coeff, const variable<RingT> &var);

    template<typename RingT>
    linear_term<RingT> operator*(const RingT &ring_coeff, const variable<RingT> &var);

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const variable<RingT> &var);

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &ring_coeff, const variable<RingT> &var);

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const variable<RingT> &var);

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &ring_coeff, const variable<RingT> &var);


/****************************** Linear term **********************************/

/**
 * A linear term represents a formal expression of the form "coeff * x_{index}".
 */
    template<typename RingT>
    class linear_term {
    public:
        var_index_t index{};
        RingT coeff;

        linear_term() = default;

        linear_term(const linear_term<RingT> &other) = default;

        linear_term(const variable<RingT> &var);

        linear_term(const variable<RingT> &var, const integer_coeff_t int_coeff);

        linear_term(const variable<RingT> &var, const RingT &ring_coeff);

        linear_term<RingT> operator*(const integer_coeff_t int_coeff) const;

        linear_term<RingT> operator*(const RingT &ring_coeff) const;

        linear_combination<RingT> operator+(const linear_combination<RingT> &other) const;

        linear_combination<RingT> operator-(const linear_combination<RingT> &other) const;

        linear_term<RingT> operator-() const;

        bool operator==(const linear_term<RingT> &other) const;
    };

    template<typename RingT>
    linear_term<RingT> operator*(const integer_coeff_t int_coeff, const linear_term<RingT> &lt);

    template<typename RingT>
    linear_term<RingT> operator*(const RingT &ring_coeff, const linear_term<RingT> &lt);

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const linear_term<RingT> &lt);

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &ring_coeff, const linear_term<RingT> &lt);

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const linear_term<RingT> &lt);

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &ring_coeff, const linear_term<RingT> &lt);


/***************************** Linear combination ****************************/

    template<typename RingT>
    class linear_combination;

    template<typename RingT>
    std::ostream &operator<<(std::ostream &out, const linear_combination<RingT> &lc);

    template<typename RingT>
    std::istream &operator>>(std::istream &in, linear_combination<RingT> &lc);

/**
 * A linear combination represents a formal expression of the form "sum_i coeff_i * x_{index_i}".
 */
    template<typename RingT>
    class linear_combination {
    public:

        std::vector<linear_term<RingT> > terms;

        linear_combination() = default;

        linear_combination(const linear_combination<RingT> &other) = default;

        linear_combination(const integer_coeff_t int_coeff);

        linear_combination(const RingT &ring_coeff);

        linear_combination(const variable<RingT> &var);

        linear_combination(const linear_term<RingT> &lt);

        linear_combination(const std::vector<linear_term<RingT> > &all_terms);

        /* for supporting range-based for loops over linear_combination */
        typename std::vector<linear_term<RingT> >::const_iterator begin() const;

        typename std::vector<linear_term<RingT> >::const_iterator end() const;

        void add_term(const variable<RingT> &var);

        void add_term(const variable<RingT> &var, const integer_coeff_t int_coeff);

        void add_term(const variable<RingT> &var, const RingT &ring_coeff);

        void add_term(const linear_term<RingT> &lt);

        RingT evaluate(const std::vector<RingT> &assignment) const;

        linear_combination<RingT> operator*(const integer_coeff_t int_coeff) const;

        linear_combination<RingT> operator*(const RingT &ring_coeff) const;

        linear_combination<RingT> operator+(const linear_combination<RingT> &other) const;

        linear_combination<RingT> operator-(const linear_combination<RingT> &other) const;

        linear_combination<RingT> operator-() const;

        bool operator==(const linear_combination<RingT> &other) const;

        [[nodiscard]] bool is_valid(size_t num_variables) const;

        void print(const std::map<size_t, std::string> &variable_annotations = std::map<size_t, std::string>()) const;

        void print_with_assignment(const std::vector<RingT> &full_assignment,
                                   const std::map<size_t, std::string> &variable_annotations = std::map<size_t, std::string>()) const;

        friend std::ostream &operator
        <<<RingT>(
        std::ostream &out,
        const linear_combination<RingT> &lc
        );

        friend std::istream &operator>><RingT>(std::istream &in, linear_combination<RingT> &lc);
    };

    template<typename RingT>
    linear_combination<RingT> operator*(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc);

    template<typename RingT>
    linear_combination<RingT> operator*(const RingT &ring_coeff, const linear_combination<RingT> &lc);

    template<typename RingT>
    linear_combination<RingT> operator+(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc);

    template<typename RingT>
    linear_combination<RingT> operator+(const RingT &ring_coeff, const linear_combination<RingT> &lc);

    template<typename RingT>
    linear_combination<RingT> operator-(const integer_coeff_t int_coeff, const linear_combination<RingT> &lc);

    template<typename RingT>
    linear_combination<RingT> operator-(const RingT &ring_coeff, const linear_combination<RingT> &lc);

} // ringsnark

#include <ringsnark/relations/variable.tcc>

#endif // VARIABLE_HPP_
