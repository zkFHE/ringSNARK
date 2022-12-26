#include <algorithm>
#include <cstdint>

using size_t = std::size_t;

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
vector<RingT> divide(const vector<RingT> &numerator, const vector<RingT> &denominator) {
    auto numerator_size = numerator.size(); // == deg(numerator) + 1
    auto denominator_size = denominator.size();

    // Invariant: numerator == denominator * quotient + remainder
    vector<RingT> quotient(numerator_size - denominator_size); // == 0 polynomial
    vector<RingT> remainder(numerator); // == numerator


    while (remainder.size() >= denominator.size()) {
        // tmp = lead(remainder) / lead(denominator)
        RingT tmp = remainder[remainder.size() - 1] / denominator[denominator.size() - 1];

        // quotient += tmp
        for (auto &q_i: quotient) {
            q_i += tmp;
        }

        // remainder -= tmp * denominator
        for (size_t i = 0; i < denominator.size(); i++) { // remainder.size() >= denominator.size()
            remainder[i] -= tmp * denominator[i];
        }

        // remove leading zero coefficients
        // after this, remainder == 0 iff remainder.size() == 0
        size_t leading_zeros = 0;
        for (int i = remainder.size() - 1; i >= 0; i--) {
            if (!remainder[i].is_zero()) { leading_zeros++; }
            else { break; }
        }
        remainder.resize(remainder.size() - leading_zeros);
    }
    return quotient;
}