#include <cassert>
#include <iostream>

#include "poly_arith.h"
#include "ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp"
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace polytools;

// #define N 8192
#define N 16384
#define LOG_T 30

/* // TODO: does not work, ciphertext modulus is not big enough for secure noise
flooding #define N 4096 #define LOG_T 16
 */

#define SEC_PARAM 128
// Most recent and smallest estimate from "Securing Approximate Homomorphic
// Encryption using Differential Privacy"
// [LMSS21](https://eprint.iacr.org/2022/816) for lambda = 128 and s = 64
#define NOISE_BITS 45

#define NUM_REPEATS 100

int main() {
  EncryptionParameters parms(scheme_type::bgv);
  size_t poly_modulus_degree = N;
  parms.set_poly_modulus_degree(poly_modulus_degree);

  // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  if (N == 8192 || N == 16384) {
    parms.set_coeff_modulus(
        CoeffModulus::Create(poly_modulus_degree, {59, 60, 60}));
  }
  parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, LOG_T));
  SEALContext context(parms);

  cout << "[PARAM] Parameter validation (success): "
       << context.parameter_error_message() << endl;
  auto qualifiers = context.first_context_data()->qualifiers();
  cout << "[PARAM] Batching enabled: " << boolalpha << qualifiers.using_batching
       << endl;
  cout << "[PARAM] poly_modulus_degree N=" << parms.poly_modulus_degree()
       << endl;
  cout << "[PARAM] plain_modulus log(t)=" << parms.plain_modulus().bit_count()
       << " bits" << endl;
  size_t coeff_modulus_bit_count = 0;
  cout << "[PARAM] coeff_modulus log(q)=";
  for (auto q_i : parms.coeff_modulus()) {
    coeff_modulus_bit_count += q_i.bit_count();
    if (q_i != parms.coeff_modulus()[0]) {
      cout << " + ";
    }
    cout << q_i.bit_count();
  }
  cout << " = " << coeff_modulus_bit_count << " bits" << endl;
  cout << "[PARAM] coeff_modulus.size=" << parms.coeff_modulus().size() << endl;

  cout << "[PARAM] plain_modulus t=" << parms.plain_modulus().value() << endl;
  cout << "[PARAM] coeff_modulus q=";
  for (auto q_i : parms.coeff_modulus()) {
    if (q_i != parms.coeff_modulus()[0]) {
      cout << " * ";
    }
    cout << q_i.value();
  }
  cout << endl;

  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);

  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);

  BatchEncoder batch_encoder(context);
  size_t slot_count = batch_encoder.slot_count();
  size_t row_size = slot_count / 2;

  auto start = chrono::high_resolution_clock::now();
  vector<Ciphertext> zeros_encrypted(SEC_PARAM);
  vector<uint64_t> zeros(slot_count, 0ULL);
  Plaintext zero_plain;
  batch_encoder.encode(zeros, zero_plain);
  for (int i = 0; i < SEC_PARAM; i++) {
    encryptor.encrypt(zero_plain, zeros_encrypted[i]);
  }
  auto end = chrono::high_resolution_clock::now();
  cout << "[TIME][CLIENT] One-time setup\t"
       << chrono::duration_cast<chrono::microseconds>(end - start).count()
       << " us" << endl;

  vector<uint64_t> pod_matrix(slot_count);
  for (size_t i = 0; i < row_size; i++) {
    pod_matrix[i] = (uint64_t)i;
    pod_matrix[row_size + i] = (uint64_t)3 * i;
  }

  start = chrono::high_resolution_clock::now();
  Plaintext x_plain;
  batch_encoder.encode(pod_matrix, x_plain);

  Ciphertext x_encrypted;
  encryptor.encrypt(x_plain, x_encrypted);
  end = chrono::high_resolution_clock::now();

  size_t time_AxR = 0;
  size_t time_RpR = 0;
  size_t time_RxR = 0;
  size_t time_enc = 0;
  size_t time_dec = 0;

  SealPoly x1(context, x_encrypted, 0);
  SealPoly x2(context, x_encrypted, 1);

  auto tables = context.get_context_data(x1.get_parms_id())->small_ntt_tables();

  // Setup for encoding space
  parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  if (N == 8192 || N == 16384) {
    parms.set_coeff_modulus(
        CoeffModulus::Create(poly_modulus_degree, {59, 60, 60}));
  }
  EncryptionParameters parms2(scheme_type::bgv);
  parms2.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 60));
  SEALContext context2(parms);
  KeyGenerator keygen2(context2);
  SecretKey secret_key2 = keygen2.secret_key();
  PublicKey public_key2;
  keygen.create_public_key(public_key2);

  Encryptor encryptor2(context2, public_key2);
  Evaluator evaluator2(context2);
  Decryptor decryptor2(context2, secret_key2);

  BatchEncoder batch_encoder2(context2);

  if (x1.is_ntt_form()) {  // Ensure x1 is in non-NTT form
    x1.intt_inplace(tables);
  }

  size_t time_NTT = 0;
  for (int i = 0; i < NUM_REPEATS; i++) {
    start = chrono::high_resolution_clock::now();
    x1.ntt_inplace(tables);
    end = chrono::high_resolution_clock::now();
    time_NTT +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();

    x1.intt_inplace(tables);
  }

  // Ensure both x1 and x2 are in NTT form
  x1.ntt_inplace(tables);
  if (!x2.is_ntt_form()) {
    x2.ntt_inplace(tables);
  }

  for (int i = 0; i < NUM_REPEATS; i++) {
    start = chrono::high_resolution_clock::now();
    x1.multiply_scalar_inplace(3);
    end = chrono::high_resolution_clock::now();
    time_AxR +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();

    start = chrono::high_resolution_clock::now();
    x1.add_inplace(x2);
    end = chrono::high_resolution_clock::now();
    time_RpR +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();

    start = chrono::high_resolution_clock::now();
    x1.multiply_inplace(x2);
    end = chrono::high_resolution_clock::now();
    time_RxR +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();

    Plaintext ptxt;
    Ciphertext ctxt;
    start = chrono::high_resolution_clock::now();
    batch_encoder2.encode(pod_matrix, ptxt);
    encryptor2.encrypt(ptxt, ctxt);
    end = chrono::high_resolution_clock::now();
    time_enc +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();

    start = chrono::high_resolution_clock::now();
    decryptor2.decrypt(ctxt, ptxt);
    batch_encoder2.decode(ptxt, pod_matrix);
    end = chrono::high_resolution_clock::now();
    time_dec +=
        chrono::duration_cast<chrono::microseconds>(end - start).count();
  }

  cout << "[TIME] NTT:   " << float(time_NTT) / NUM_REPEATS << " us" << endl;
  cout << "[TIME] A x R: " << float(time_AxR) / NUM_REPEATS << " us" << endl;
  cout << "[TIME] R + R: " << float(time_RpR) / NUM_REPEATS << " us" << endl;
  cout << "[TIME] R x R: " << float(time_RxR) / NUM_REPEATS << " us" << endl;
  cout << "[TIME] 1 Enc: " << float(time_enc) / NUM_REPEATS << " us" << endl;
  cout << "[TIME] 1 Dec: " << float(time_dec) / NUM_REPEATS << " us" << endl;

  unsigned long long size = 9 * x_encrypted.size() *
                            parms.coeff_modulus().size() *
                            parms.poly_modulus_degree() * 8;
  cout << "[SPACE] Proof size\t" << size << " B" << endl;
  cout << endl;
}
