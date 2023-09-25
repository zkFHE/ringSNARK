#include <benchmark/benchmark.h>

#include <ringsnark/gadgetlib/protoboard.hpp>
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include "ringsnark/seal/seal_util.hpp"
#include "seal/seal.h"

typedef ringsnark::seal::RingElem R;
typedef ringsnark::seal::EncodingElem E;
using namespace std;
using namespace seal;
using namespace ringsnark;

#define USE_MODSWITCH_IN_INNER_PRODUCT

// Set such that the noise in the outer computation does not overflow
const size_t poly_modulus_degree = 2048;

// Set such that the noise in the inner computation (within the proof) does not overflow
const size_t blowup_poly_modulus_degree = 8;
const size_t inner_poly_modulus_degree = blowup_poly_modulus_degree * poly_modulus_degree;

const size_t logT = 16;
const size_t num_features = 256;

SEALContext init_context() {
  EncryptionParameters params(scheme_type::bgv);

  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_coeff_modulus(default_double_batching_modulus(
      poly_modulus_degree, inner_poly_modulus_degree));
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, logT));
  SEALContext context(params);

  return context;
}

std::pair<vector<Ciphertext>, vector<Ciphertext>> encrypt() {
  auto context = init_context();
  BatchEncoder encoder(context);
  KeyGenerator keygen(context);
  Encryptor encryptor(context, keygen.secret_key());

  vector<uint64_t> vs(poly_modulus_degree);
  Plaintext ptxt;
  vector<Ciphertext> ctxt1(num_features), ctxt2(num_features);

  for (size_t i = 0; i < num_features; i++) {
    vs = {i};
    encoder.encode(vs, ptxt);
    encryptor.encrypt_symmetric(ptxt, ctxt1[i]);
    assert(ctxt1[i].is_ntt_form());

    vs = {2*i+1 % (1<<logT)};
    encoder.encode(vs, ptxt);
    encryptor.encrypt_symmetric(ptxt, ctxt2[i]);
    assert(ctxt2[i].is_ntt_form());
  }
  return std::make_pair(ctxt1, ctxt2);
}

ringsnark::protoboard<R> init(SEALContext& context) {
  try {
    R::set_context(context);
    E::set_context(inner_poly_modulus_degree);
  } catch (std::invalid_argument &e) {
  }

  const size_t num_input_variables = 2 * num_features + 5;
  ringsnark::protoboard<R> pb;

  vector<pb_variable_array<R>> in1(num_features), in2(num_features);
  for (size_t i = 0; i < num_features; i++) {
    in1[i].allocate(pb, 2, "in1" + to_string(i));
    in2[i].allocate(pb, 2, "in2" + to_string(i));
  }
  pb_variable_array<R> out;
  out.allocate(pb, 5, "out");
  pb.set_input_sizes(num_input_variables);

  vector<linear_combination<R>> sum(3);
  pb_variable_array<R> prods_00, prods_01, prods_10, prods_11;
  prods_00.allocate(pb, num_features, "prods_00");
  prods_01.allocate(pb, num_features, "prods_01");
  prods_10.allocate(pb, num_features, "prods_10");
  prods_11.allocate(pb, num_features, "prods_11");

  // multiply and aggregate
  for (size_t i = 0; i < num_features; i++) {
    pb.add_r1cs_constraint(
        r1cs_constraint<R>(in1[i][0], in2[i][0], prods_00[i]));
    sum[0] = sum[0] + prods_00[i];

    pb.add_r1cs_constraint(
        r1cs_constraint<R>(in1[i][0], in2[i][1], prods_01[i]));
    pb.add_r1cs_constraint(
        r1cs_constraint<R>(in1[i][1], in2[i][0], prods_10[i]));
    sum[1] = sum[1] + prods_01[i] + prods_10[i];

    pb.add_r1cs_constraint(
        r1cs_constraint<R>(in1[i][1], in2[i][1], prods_11[i]));
    sum[2] = sum[2] + prods_11[i];
  }

  // approx sigmoid of degree 2 (i.e., do a squaring)
  vector<pb_linear_combination<R>> sq(5);

  pb_variable<R> s02, s11;
  s02.allocate(pb, "s02");
  s11.allocate(pb, "s11");

  pb.add_r1cs_constraint(r1cs_constraint<R>(sum[0], sum[0], out[0]));

  pb.add_r1cs_constraint(r1cs_constraint<R>(2 * sum[0], sum[1], out[1]));

  pb.add_r1cs_constraint(r1cs_constraint<R>(sum[0], sum[2], s02));
  pb.add_r1cs_constraint(r1cs_constraint<R>(sum[1], sum[1], s11));
  pb.add_r1cs_constraint(r1cs_constraint<R>(1, 2 * s02 + s11, out[2]));

  pb.add_r1cs_constraint(r1cs_constraint<R>(sum[1], sum[2], out[3]));

  pb.add_r1cs_constraint(r1cs_constraint<R>(sum[2], sum[2], out[4]));

  cout << "#inputs:      " << pb.num_inputs() << endl;
  cout << "#variables:   " << pb.num_variables() << endl;
  cout << "#constraints: " << pb.num_constraints() << endl;
  cout << endl;

  // Set values
  {
    BatchEncoder encoder(context);
    KeyGenerator keygen(context);
    Encryptor encryptor(context, keygen.secret_key());
    Evaluator evaluator(context);

    vector<uint64_t> vs(poly_modulus_degree);
    Plaintext ptxt;
    //    vector<Ciphertext> ctxt1(num_features), ctxt2(num_features);

    Ciphertext sum_ctxt;
    const auto &[ctxt1, ctxt2] = encrypt();
    assert(ctxt1.size() == num_features);
    assert(ctxt2.size() == num_features);
    for (size_t i = 0; i < num_features; i++) {
      // Set inputs
      assert(ctxt1[i].is_ntt_form());
      const auto in1_p0 = R(polytools::SealPoly(context, ctxt1[i], 0));
      const auto in1_p1 = R(polytools::SealPoly(context, ctxt1[i], 1));
      pb.val(in1[i][0]) = in1_p0;
      pb.val(in1[i][1]) = in1_p1;

      assert(ctxt2[i].is_ntt_form());
      const auto in2_p0 = R(polytools::SealPoly(context, ctxt2[i], 0));
      const auto in2_p1 = R(polytools::SealPoly(context, ctxt2[i], 1));
      pb.val(in2[i][0]) = in2_p0;
      pb.val(in2[i][1]) = in2_p1;

      // Set products
      // Because of the R1CS arithmetization, we need access to the 0-1 and 1-0
      // products separately
      const auto p00 = in1_p0 * in2_p0;
      const auto p01 = in1_p0 * in2_p1;
      const auto p10 = in1_p1 * in2_p0;
      const auto p11 = in1_p1 * in2_p1;

      pb.val(prods_00[i]) = p00;
      pb.val(prods_01[i]) = p01;
      pb.val(prods_10[i]) = p10;
      pb.val(prods_11[i]) = p11;

      // Aggregate
      Ciphertext prod;
      evaluator.multiply(ctxt1[i], ctxt2[i], prod);
      if (i == 0) {
        sum_ctxt = prod;
      } else {
        evaluator.add_inplace(sum_ctxt, prod);
      }
    }

    // Apply sigmoid approximation of degree 2
    const auto s0 = R(polytools::SealPoly(context, sum_ctxt, 0));
    const auto s1 = R(polytools::SealPoly(context, sum_ctxt, 1));
    const auto s2 = R(polytools::SealPoly(context, sum_ctxt, 2));

    const auto t0 = s0 * s0;
    pb.val(out[0]) = t0;

    const auto t1 = 2 * s0 * s1;
    pb.val(out[1]) = t1;

    const auto s_02 = s0 * s2;
    pb.val(s02) = s_02;
    const auto s_11 = s1 * s1;
    pb.val(s11) = s_11;
    const auto t2 = 2 * s_02 + s_11;
    pb.val(out[2]) = t2;

    const auto t3 = s1 * s2;
    pb.val(out[3]) = t3;

    const auto t4 = s2 * s2;
    pb.val(out[4]) = t4;
  }

  cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl;
  cout << endl;

  return pb;
}

Ciphertext eval(const Evaluator &evaluator, const vector<Ciphertext> &ctxt1,
                const vector<Ciphertext> &ctxt2) {
  Ciphertext out;
  assert(ctxt1.size() == num_features);
  assert(ctxt2.size() == num_features);
  Ciphertext sum;
  Ciphertext prod;
  for (size_t i = 0; i < num_features; i++) {
    evaluator.multiply(ctxt1[i], ctxt2[i], prod);
    if (i == 0) {
      sum = prod;
    } else {
      evaluator.add_inplace(sum, prod);
    }
  }
  evaluator.square_inplace(sum);
  return sum;
}

ringsnark::protoboard<R> *global_protoboard;
ringsnark::rinocchio::keypair<R, E> *global_keypair;
ringsnark::rinocchio::proof<R, E> *global_proof;

static void BM_LogRegInf_Setup(benchmark::State &state) {
  ringsnark::protoboard<R> pb;
  if (!global_protoboard) {
    auto context = init_context();
    print_params(context.first_context_data()->parms());
    global_protoboard = new ringsnark::protoboard<R>(init(context));
  }
  pb = *global_protoboard;

  for (auto _ : state) {
    auto kp = ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
    if (!global_keypair) {
      global_keypair = new ringsnark::rinocchio::keypair<R, E>(kp);
    }
  }
}

static void BM_LogRegInf_Prove(benchmark::State &state) {
  const auto pb = *global_protoboard;
  const auto kp = *global_keypair;
  for (auto _ : state) {
    auto pi = ringsnark::rinocchio::prover(kp.pk, pb.primary_input(),
                                           pb.auxiliary_input());
    if (!global_proof) {
      global_proof = new ringsnark::rinocchio::proof<R, E>(pi);
    }
  }
}

static void BM_LogRegInf_Verify(benchmark::State &state) {
  const auto pb = *global_protoboard;
  const auto kp = *global_keypair;
  const auto pi = *global_proof;
  for (auto _ : state) {
    bool verif = ringsnark::rinocchio::verifier(kp.vk, pb.primary_input(), pi);
    benchmark::DoNotOptimize(verif);
  }
}

static void BM_LogRegInf_Keygen(benchmark::State &state) {
  auto context = init_context();
  for (auto _ : state) {
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    benchmark::DoNotOptimize(secret_key);
  }
}

static void BM_LogRegInf_Eval(benchmark::State &state) {
  auto context = init_context();
  const auto &[ctxt1, ctxt2] = encrypt();
  Evaluator evaluator(context);
  Ciphertext out;
  for (auto _ : state) {
    out = eval(evaluator, ctxt1, ctxt2);
    benchmark::DoNotOptimize(out);
  }
}

static void BM_LogRegInf_EncDec(benchmark::State &state) {
  auto context = init_context();
  const auto &[c1, c2] = encrypt();
  Evaluator evaluator(context);
  KeyGenerator keygen(context);
  Decryptor decryptor(context, keygen.secret_key());
  const auto out = eval(evaluator, c1, c2);
  Plaintext ptxt;
  for (auto _ : state) {
    auto [ctxt1, ctxt2] = encrypt();
    benchmark::DoNotOptimize(ctxt1);
    benchmark::DoNotOptimize(ctxt2);

    decryptor.decrypt(out, ptxt);
    benchmark::DoNotOptimize(ptxt);
  }
}

BENCHMARK(BM_LogRegInf_Setup);
BENCHMARK(BM_LogRegInf_Prove);
BENCHMARK(BM_LogRegInf_Verify);

BENCHMARK(BM_LogRegInf_Keygen);
BENCHMARK(BM_LogRegInf_Eval);
BENCHMARK(BM_LogRegInf_EncDec);

BENCHMARK_MAIN();