#include <cmath>
#include <iostream>

#include "seal/seal.h"
#include "ring.h"
#include "ring_snark.h"
#include "poly_arith.h"

using namespace std;
using namespace seal;
using namespace rinocchio;
using E = SealPoly;
using R = SealPoly;
using A = uint64_t;

// Circuit is x6=(x1+x2)*x5; x5=(x3*x4)
bool prove_circuit(SEALContext &context, EncryptionParameters &params, vector<SealPoly> &values_io,
                   vector<SealPoly> &values_mid) {
    assert(all_of(values_io.begin(), values_io.end(), [](const E &e) { return e.is_ntt_form(); }));
    assert(all_of(values_mid.begin(), values_mid.end(), [](const E &e) { return e.is_ntt_form(); }));

    // Parameters
    // Randomly pick r5, r6 <- A
    uint64_t q1 = params.coeff_modulus()[0].value();
    uint64_t r5, r6;
    r5 = 2;//rand() % (p1.value() - 1) + 1;
    do {
        r6 = 3; //rand() % (p1.value() - 1) + 1;
    } while (r5 == r6);
    uint64_t inv = 1; // (r6-r5)^-1 = 1, since (3-2) * 1 == 1 mod q
    uint64_t min_inv = (q1 - inv) % q1;
    uint64_t r5_inv = (r5 * inv) % q1;
    uint64_t min_r5_inv = (q1 - r5_inv) % q1;
    uint64_t r5_r6 = (r5 * r6) % q1;

    R c1 = values_io[0];
    R c2 = values_io[1];
    R c3 = values_io[2];
    R c4 = values_io[3];
    R c5 = values_mid[0];
    R c6 = values_io[4];

    auto start = chrono::high_resolution_clock::now();
    vector<vector<A>> v{
            vector<A>{min_r5_inv, inv},
            vector<A>{min_r5_inv, inv},
            vector<A>{(1 + r5_inv) % q1, min_inv},
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{0, 0}
    };

    vector<vector<A>> w{
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{(1 + r5_inv) % q1, min_inv},
            vector<A>{min_r5_inv, inv},
            vector<A>{0, 0}
    };

    vector<vector<A>> y{
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{0, 0},
            vector<A>{(1 + r5_inv) % q1, min_inv},
            vector<A>{min_r5_inv, inv}
    };

    uint64_t min_r5_min_r6 = (q1 - r5) % q1;
    min_r5_min_r6 = (q1 + min_r5_min_r6 - r6) % q1;
    vector<A> t{r5_r6, min_r5_min_r6, 1};

    //// h: [4 * ((r6-r5)^-1)^2 * (c1+c2-c3) * (c5-c4)]
    // h: ((r6-r5)^-1)^2 * (c1+c2-c3) * (c1+c2-c4)
    R h0(c1);
    h0.add_inplace(c2);
    R tmp(h0);
    h0.subtract_inplace(c3);
    tmp.subtract_inplace(c4);

    h0.multiply_inplace(tmp);

    uint64_t scalar = (inv * inv) % q1;
    h0.multiply_scalar_inplace(scalar);

    Poly<R, A> h{vector<R>{h0}};

    vector<size_t> indices_io{0, 1, 2, 3, 5};
    vector<size_t> indices_mid{4};
    vector<R> values{c1, c2, c3, c4, c5, c6};
    SnarkParameters<R, A> snarkParameters(context, values, indices_io, indices_mid, v, w, y, t, h);

    Prover<E, R, A> prover;
    Verifier<E, R, A> verifier(snarkParameters, params);

    // TODO
    const R &alpha(c1);

    A s = 3;
    const R &beta(c2);

    R r_v(c3);
    const R &r_w(c4);
    //r_v.multiply_scalar_inplace(0);
    R r_y(r_v);
    r_y.multiply_inplace(r_w);

    auto end = chrono::high_resolution_clock::now();
    cout << "=======================================" << endl;
    cout << "QRP setup\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl;

    start = chrono::high_resolution_clock::now();
    auto vk = verifier.generate_vk(snarkParameters, s, alpha, beta, r_v, r_w, r_y);
    end = chrono::high_resolution_clock::now();
    cout << "Setup (CRS, vk)\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl;


    auto crs = vk.crs;

    start = chrono::high_resolution_clock::now();
    auto proof = prover.prove(snarkParameters, crs);
    end = chrono::high_resolution_clock::now();
    cout << "Prove\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl;

    start = chrono::high_resolution_clock::now();
    bool verified = verifier.verify(vk, snarkParameters, proof);
    end = chrono::high_resolution_clock::now();
    cout << "Verify\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl << flush;
    return verified;
}

int main() {
    EncryptionParameters params(scheme_type::bfv);
    auto poly_modulus_degree = (size_t) pow(2, 12);
    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(params);

    cout << "Parameter validation (success): "
         << context.parameter_error_message() << endl;
    auto qualifiers = context.first_context_data()->qualifiers();
    cout << "Batching enabled: " << boolalpha << qualifiers.using_batching
         << endl;

    KeyGenerator keygen(context);
    const SecretKey &secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    BatchEncoder batch_encoder(context);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    size_t slot_count = batch_encoder.slot_count();
    size_t row_size = slot_count / 2;

    vector<uint64_t> pod_matrix(slot_count, 0ULL);
    for (size_t i = 0; i < row_size; i++) {
        pod_matrix[i] = (uint64_t) i;
        pod_matrix[row_size + i] = (uint64_t) 3 * i;
    }

    Plaintext x_plain;
    batch_encoder.encode(pod_matrix, x_plain);

    for (size_t i = 0; i < row_size; i++) {
        pod_matrix[i] = (uint64_t) i + 1;
        pod_matrix[row_size + i] = (uint64_t) 2 * i;
    }

    Plaintext y_plain;
    batch_encoder.encode(pod_matrix, y_plain);

    Ciphertext x_enc, y_enc;
    encryptor.encrypt(x_plain, x_enc);
    encryptor.encrypt(y_plain, y_enc);


    auto start = chrono::high_resolution_clock::now();

    SealPoly x1(context, x_enc, 0);
    SealPoly x2(context, x_enc, 1);
    SealPoly x3(context, y_enc, 0);
    SealPoly x4(context, y_enc, 1);

    auto tables = context.get_context_data(x1.get_parms_id())->small_ntt_tables();

    x1.ntt_inplace(tables);
    x2.ntt_inplace(tables);
    x3.ntt_inplace(tables);
    x4.ntt_inplace(tables);
    auto end = chrono::high_resolution_clock::now();
    cout << "=======================================" << endl;
    cout << "Eval: 4 NTT\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl;

    start = chrono::high_resolution_clock::now();

    SealPoly x5(x3);
    x5.multiply_inplace(x4);

    SealPoly tmp(x1);
    tmp.add_inplace(x2);

    SealPoly x6(x5);
    x6.multiply_inplace(tmp);

    end = chrono::high_resolution_clock::now();
    cout << "=======================================" << endl;
    cout << "Eval{x6=(x1+x2)*x5; x5=x3*x4}\t" << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us"
         << endl;

    vector<SealPoly> values_io{x1, x2, x3, x4, x6};
    vector<SealPoly> values_mid{x5};

    bool verified = prove_circuit(context, params, values_io, values_mid);
    cout << "Verification: " << verified << endl;


    /**/
//    Plaintext x_ptxt;
//    decryptor.decrypt(x_enc, x_ptxt);
//    vector<uint64_t> x_vec;
//    batch_encoder.decode(x_ptxt, x_vec);
//    cout << x_vec[0] << endl;
//
//    Plaintext zero;
//    vector<uint64_t> zero_vec(params.poly_modulus_degree());
//    batch_encoder.encode(zero_vec, zero);
//    evaluator.multiply_plain_inplace(x_enc, zero);
//    decryptor.decrypt(x_enc, x_ptxt);
//    vector<uint64_t> x_vec2;
//    batch_encoder.decode(x_ptxt, x_vec2);
//    cout << x_vec2[0] << endl;

    return verified;
}