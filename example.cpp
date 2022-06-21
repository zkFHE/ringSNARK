#include <cmath>
#include <iostream>

#include "seal/seal.h"
#include "ring.h"
#include "poly_arith.h"

using namespace std;
using namespace seal;
using namespace rinocchio;

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

    vector <uint64_t> pod_matrix(slot_count, 0ULL);
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


    Ciphertext res;
    auto start = chrono::high_resolution_clock::now();
    evaluator.multiply(x_enc, y_enc, res);
    auto end = chrono::high_resolution_clock::now();
    cout << "=======================================" << endl;
    cout << "Mult: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " us" << endl;


    Plaintext decrypted_result;
    decryptor.decrypt(res, decrypted_result);
    vector <uint64_t> res_vec;
    batch_encoder.decode(decrypted_result, res_vec);

    return 0;
}