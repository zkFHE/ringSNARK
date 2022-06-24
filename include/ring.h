#ifndef RING_H
#define RING_H

#include <iostream>
#include <vector>
#include "seal/seal.h"
#include "poly_arith.h"

using namespace polytools;
using namespace std;

namespace rinocchio {
    // R is the underlying ring, A is the set of exceptional elements of R
    template<typename R, typename A>
    class Poly {
    public:
        vector<R> coeffs;

        static Poly<R, A> zero() {
            return Poly<R, A>();
        }

        static Poly<R, A> one() {
            return Poly<R, A>();
        }

        static Poly<R, A> random_element() {
            return Poly<R, A>();
        }

        Poly<R, A> operator-() const {
            return Poly<R, A>();
        }

        Poly<R, A> operator+(Poly<R, A> const &other) {
            return Poly<R, A>();
        }

        Poly<R, A> &operator+=(Poly<R, A> const &other) {
            return Poly<R, A>();
        }

        Poly<R, A> operator-(Poly<R, A> const &other) {
            return Poly<R, A>();
        }

        Poly<R, A> operator*(Poly<R, A> const &other) {
            return Poly<R, A>();
        }

        Poly<R, A> operator*(A const &other) {
            return Poly<R, A>();
        }

        R eval_at(A x) {
            R res = R(coeffs[coeffs.size() - 1]);
            if (coeffs.size() > 1) { // deg >= 1
                for (int i = coeffs.size() - 2; i >= 0; i--) {
                    res.multiply_scalar_inplace(x);
                    res.add_inplace(coeffs[i]);
                }
            }
            return res;
        }
    };
}
#endif