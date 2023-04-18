#include "polynomials.hpp"
#include <unordered_set>
#include <stdexcept>

namespace ringsnark {
    template<typename RingT>
    evaluation_domain<RingT>::evaluation_domain(const size_t m) : m(m), values(vector<RingT>(m)) {
        for (uint64_t i = 0; i < m; i++) {
            values[i] = RingT(i); // TODO: assert i is in exceptional set
        }
    }

    template<typename RingT>
    RingT evaluation_domain<RingT>::get_domain_element(const size_t idx) const {
        return values[idx];
    }

    template<typename RingT>
    vector<RingT> evaluation_domain<RingT>::evaluate_all_lagrange_polynomials(const RingT &t) const {
        if (std::find(values.begin(), values.end(), t) != values.end()) {
            throw std::invalid_argument("t cannot be one of the values in the domain");
        }
        // TODO: Optimize
        // lagrange[j] = \prod_{i=0,i!=j}^m (t-x_i) / (x_j-x_i)
        vector<RingT> lagrange(m);
        for (int j = 0; j < m; j++) {
            lagrange[j] = RingT::one();
            RingT denominator = RingT::one();
            for (int i = 0; i < m; i++) {
                if (i != j) {
                    lagrange[j] *= (t - values[i]);
                    denominator *= (values[j] - values[i]);
                }
            }
            lagrange[j] /= denominator;
        }
        return lagrange;
    }

    template<typename RingT>
    RingT evaluation_domain<RingT>::compute_vanishing_polynomial(const RingT &t) const {
        RingT res = t - values[0];
        for (size_t i = 1; i < m; i++) {
            res *= (t - values[i]);
        }
        return res;
    }

    template<typename RingT>
    vector<RingT> evaluation_domain<RingT>::vanishing_polynomial() const {
        polynomial<RingT> Z_poly(vector<RingT>{-values[0], RingT::one()});
        for (size_t i = 1; i < values.size(); i++) {
            Z_poly *= polynomial<RingT>(vector<RingT>{-values[i], RingT::one()});
        }
        return Z_poly.data();
    }

    template<typename RingT>
    void evaluation_domain<RingT>::add_poly_Z(const RingT &coeff, vector<RingT> &H) const {
        vector<RingT> Z = vanishing_polynomial();

        for (size_t i = 0; i < std::min(H.size(), Z.size()); i++) {
            H[i] += coeff * Z[i];
        }
        if (H.size() < Z.size()) {
            H.resize(Z.size());
            for (size_t i = H.size(); i < Z.size(); i++) {
                H[i] = coeff * Z[i];
            }
        }

    }

    // TODO: not sure if we are technically dividing on coset, we might be using a misleading method name at the moment.
    template<typename RingT>
    void evaluation_domain<RingT>::divide_by_Z_on_coset(vector<RingT> &P) const {
        // implementation also normalizes P, i.e., strips zero higher coefficients
        P = divide(P, vanishing_polynomial());
    }

    template<typename RingT>
    shared_ptr<evaluation_domain<RingT>> get_evaluation_domain(const size_t min_size) {
        shared_ptr<evaluation_domain<RingT>>
                shared_domain(new evaluation_domain<RingT>(min_size));
        return shared_domain;
    }
}