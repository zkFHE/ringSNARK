#include <algorithm>
#include <cstdint>
#include <boost/math/tools/polynomial.hpp>

using size_t = std::size_t;
template<typename RingT>
using polynomial = boost::math::tools::polynomial<RingT>;

template<typename RingT>
vector<RingT> interpolate(const vector<RingT> &x, const vector<RingT> &y) {
    assert(x.size() == y.size());
    int n = x.size();
    int k, j, i;
    RingT phi, ff, b;

    vector<RingT> coeffs(n);
    vector<RingT> s(n);
    for (i = 0; i < n; i++) {
        s[i] = RingT::zero();
        coeffs[i] = RingT::zero();
    }
    s[n - 1] = -x[0];

    for (i = 1; i < n; i++) {
        for (j = n - i - 1; j < n - 1; j++) {
            s[j] -= x[i] * s[j + 1];
        }
        s[n - 1] -= x[i];
    }
    for (j = 0; j < n; j++) {
        phi = RingT(n);
        for (k = n - 1; k > 0; k--) {
            phi = RingT(k) * s[k] + x[j] * phi; // TODO: use *= and break down
        }
        ff = y[j] / phi;
        b = RingT::one();
        for (k = n - 1; k >= 0; k--) {
            coeffs[k] += b * ff;
            b = s[k] + x[j] * b;
        }
    }
    return coeffs;
}

template<typename RingT>
RingT eval(const vector<RingT> &coeffs, const RingT &x) {
    RingT res(coeffs.size() - 1);
    for (int i = coeffs.size() - 2; i >= 0; i--) {
        res *= x;
        res += coeffs[i];
    }
    return res;
}

template<typename RingT>
inline bool is_zero(const vector<RingT> &coeffs) {
    return std::all_of(coeffs.begin(), coeffs.end(), [](const RingT &c) { return c.is_zero(); });
}

template<typename RingT>
vector<RingT> multiply(const vector<RingT> &x, const vector<RingT> &y) {
    polynomial<RingT> x_poly(x), y_poly(y);
    x_poly *= y_poly;
    return x_poly.data();

}

template<typename RingT>
vector<RingT> add(const vector<RingT> &x, const vector<RingT> &y) {
    polynomial<RingT> x_poly(x), y_poly(y);
    x_poly += y_poly;
    return x_poly.data();
}

template<typename RingT>
vector<RingT> divide(const vector<RingT> &numerator, const vector<RingT> &denominator) {
    polynomial<RingT> x_poly(numerator), y_poly(denominator);
    x_poly /= y_poly;
    return x_poly.data();
}