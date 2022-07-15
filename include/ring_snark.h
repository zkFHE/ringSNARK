#ifndef RING_SNARK_H
#define RING_SNARK_H

#include <utility>
#include <set>
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
                        const vector<vector<A>> &v, const vector<vector<A>> &w, const vector<vector<A>> &y,
                        const vector<A> &t, const Poly<R, A> &h) :
                context(std::move(context)), values(values),
                indices_io(indices_io), indices_mid(indices_mid),
                v(v), w(w), y(y), t(t), h(h) {
            num_io = indices_io.size();
            num_mid = indices_mid.size();
            values_mid = vector<R>();
            set<size_t> indices;
            if (v.size() != w.size() || v.size() != y.size()) {
                throw invalid_argument("v, w, y must have the same size");
            }
            size_t max_index = v.size() - 1;

            for (size_t i_mid: indices_mid) {
                if (i_mid > max_index) {
                    throw invalid_argument("intermediate index too big");
                }
                values_mid.push_back(values[i_mid]);
                if (indices.count(i_mid)) {
                    throw invalid_argument("duplicate intermediate index " + to_string(i_mid));
                }
                indices.insert(i_mid);
            }
            values_io = vector<R>();
            for (size_t i_io: indices_io) {
                if (i_io > max_index) {
                    throw invalid_argument("input/output index too big");
                }
                values_io.push_back(values[i_io]);
                if (indices.count(i_io)) {
                    throw invalid_argument("duplicate input/output index " + to_string(i_io));
                }
                indices.insert(i_io);
            }
            if (indices.size() != v.size()) {
                throw invalid_argument("mismatched number of indices and number of polynomials");
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
        Prover() = default;

        ~Prover() = default;

        Proof<E> prove(SnarkParameters<R, A> params, CRS<E> crs);

    protected:
        E inner_prod_enc_ring(const vector<R> &poly, const vector<E> &values);

        E inner_prod_enc_exc(const vector<A> &poly, const vector<E> &values);

        virtual E multiply_inplace(E &e, const A &a) = 0;

        virtual E multiply_inplace(E &e, const R &r) = 0;

        virtual E add_inplace(E &e1, const E &e2) = 0;
    };

    template<typename E, typename R, typename A>
    class Verifier {
    public:
        SnarkParameters<R, A> snark_parameters;

        explicit Verifier(SnarkParameters<R, A> snark_parameters) : snark_parameters(snark_parameters) {}


        VerifKey<E, R, A> generate_vk(SnarkParameters<R, A> snark_parameters, A s, R alpha, R beta, R r_v,
                                      R r_w, R r_y);

        bool verify(VerifKey<E, R, A> vk, SnarkParameters<R, A> params, Proof<E> proof);

    protected:
        A eval(const vector<A> &poly, A x);

        virtual R decode(void *sk, const E &x) = 0;

        virtual E encode(void *pk, const R &x) = 0;

        virtual E encode(void *pk, const A &x) = 0;

        virtual A multiply_inplace(A &a1, const A &a2) = 0;

        virtual R multiply_inplace(R &r, const A &a) = 0;

        virtual R multiply_inplace(R &r1, const R &r2) = 0;

        virtual A add_inplace(A &a1, const A &a2) = 0;

        virtual R add_inplace(R &r1, const R &r2) = 0;

        virtual R subtract_inplace(R &r1, const R &r2) = 0;

        virtual bool is_unit(const R &x) = 0;

        virtual bool is_zero(const A &a) = 0;

        virtual bool is_zero(const R &r) = 0;

        virtual bool is_equal(const A &a1, const A &a2) = 0;

        virtual bool is_equal(const R &r1, const R &r2) = 0;
    };

    template<typename R, typename A>
    SnarkParameters<R, A> setup(uint64_t security_param);
}  // namespace rinocchio

#include "prover.tpp"
#include "verifier.tpp"

#endif