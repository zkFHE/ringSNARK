#ifndef RING_SNARK_H
#define RING_SNARK_H

#include <cmath>
#include <iostream>
#include <utility>
#include "ring.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

namespace rinocchio {
    template<typename R, typename A>
    class SnarkParameters {
    public:
        SnarkParameters(SEALContext context, const vector<R> &values,
                        const vector<size_t> &indices_io, const vector<size_t> &indices_mid,
                        const Poly<R, A> &v0, const Poly<R, A> &w0,
                        const Poly<R, A> &y0, const vector<vector<A>> &v, const vector<vector<A>> &w,
                        const vector<vector<A>> &y,
                        const vector<A> &t, const Poly<R, A> &h) :
                context(std::move(context)), values(values),
                indices_io(indices_io), indices_mid(indices_mid), v0(v0), w0(w0), y0(y0),
                v(v), w(w), y(y), t(t), h(h) {
            num_io = indices_io.size();
            num_mid = indices_mid.size();
            values_mid = vector<R>();
            for (size_t i_mid: indices_mid) {
                values_mid.push_back(values[i_mid]);
            }
            values_io = vector<R>();
            for (size_t i_io: indices_io) {
                values_io.push_back(values[i_io]);
            }
        }

        SEALContext context;
        vector<R> values;
        vector<R> values_mid;
        vector<R> values_io;
        vector<size_t> indices_io;
        vector<size_t> indices_mid;
        size_t num_io;
        size_t num_mid;
        Poly<R, A> v0, w0, y0;
        vector<vector<A>> v, w, y;
        vector<A> t;
        Poly<R, A> h;
    };

    // E is the encoding space
    template<typename E>
    struct CRS {
        vector<E> s_pows, alpha_s_pows;    // size = |I_mid|
        vector<E> beta_prod;               // size = |I_mid|
        void *pk;
    };

    template<typename E>
    struct Proof {
        E A;
        E A_prime;
        E B;
        E B_prime;
        E C;
        E C_prime;
        E D;
        E D_prime;
        E F;
    };

    template<typename E, typename R, typename A>
    struct VerifKey {
        VerifKey() = default;

        void *sk{};
        CRS<E> crs;
        A s;
        R alpha;
        R beta;
        R r_v, r_w, r_y;
    };

    template<typename E, typename R, typename A>
    class Prover {
    public:
        Proof<E> prove(SnarkParameters<R, A> params, CRS<E> crs);

        E sum(const vector<E> &values);
        E inner_prod_enc_ring(const vector<R> &poly, const vector<E> &values);
        E inner_prod_enc_exc(const vector<A> &poly, const vector<E> &values);
    };

    template<typename E, typename R, typename A>
    class Verifier {
    public:
        SnarkParameters<R, A> snark_parameters;
        EncryptionParameters encryption_parameters;

        explicit Verifier(SnarkParameters<R, A> snark_parameters, EncryptionParameters encryption_parameters)
                : snark_parameters(snark_parameters),
                  encryption_parameters(encryption_parameters) {}

        VerifKey<E, R, A> generate_vk(SnarkParameters<R, A> snark_parameters, A s, R alpha, R beta, R r_v,
                                      R r_w, R r_y);

        bool verify(VerifKey<E, R, A> vk, SnarkParameters<R, A> params, Proof<E> proof);

        A eval(const vector<A> &poly, const A x);

        R decode(void *sk, E x);

        E encode(void *pk, R x);

        E encode(void *pk, A x);
    };

    template<typename R, typename A>
    SnarkParameters<R, A> setup(uint64_t security_param);
}  // namespace rinocchio

#include "ring_snark.tpp"

#endif