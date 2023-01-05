#ifndef RINGSNARK_INTERPOLATION_H
#define RINGSNARK_INTERPOLATION_H

#include <vector>

using std::vector;

/**
 *
 * @tparam RingT
 * @param x vector of n evaluation points
 * @param x vector of n evaluation points
 * @param y vector of n evaluations y[i] = f(x[i]), for a degree-(n-1) polynomial f
 * @return coefficients `coeffs' of f, such that y[i] = sum_j coeffs[j] * x[i]^j
 */
template<typename RingT>
vector<RingT> interpolate(const vector<RingT> &x, const vector<RingT> &y);

template<typename RingT>
RingT eval(const vector<RingT> &coeffs, const RingT &x);

template<typename RingT>
vector<RingT> multiply(const vector<RingT> &x, const vector<RingT> &y);

template<typename RingT>
vector<RingT> add(const vector<RingT> &x, const vector<RingT> &y);

/**
 * Return quotient of polynomial division of `numerator' by `denominator' over the ring RingT.
 * The result is a vector of coefficients of size numerator.size() - denominator.size() (padded with zeros if the
 * quotient polynomial has a degree less than this size).
 * @tparam RingT
 * @param numerator
 * @param denominator
 * @return
 */
template<typename RingT>
vector<RingT> divide(const vector<RingT> &numerator, const vector<RingT> &denominator);

#include "polynomials.tcc"

#endif //RINGSNARK_INTERPOLATION_H
