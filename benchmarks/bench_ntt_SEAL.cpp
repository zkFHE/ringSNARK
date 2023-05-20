#include <benchmark/benchmark.h>

#include <ringsnark/gadgetlib/protoboard.hpp>
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include "seal/seal.h"

typedef ringsnark::seal::RingElem R;
typedef ringsnark::seal::EncodingElem E;

using namespace std;
using namespace seal;

ringsnark::protoboard<R> init(size_t N) {
  EncryptionParameters params(scheme_type::bgv);
  params.set_poly_modulus_degree(N);
  params.set_coeff_modulus(CoeffModulus::BFVDefault(N));
  params.set_plain_modulus(PlainModulus::Batching(N, 16));
  SEALContext context(params);

  try {
    R::set_context(context);
    E::set_context();
  } catch (std::invalid_argument &e) {
  }

  const size_t num_input_variables = N + 1;
  const size_t num_variables = N + 1;
  ringsnark::protoboard<R> pb;

  // Set public values
  ringsnark::pb_variable_array<R> vars(num_variables,
                                       ringsnark::pb_variable<R>());
  vars.allocate(pb, num_variables, "x");
  pb.set_input_sizes(num_input_variables);

  // Set constraints
  auto tables =
      context.get_context_data(context.first_parms_id())->small_ntt_tables();
  vector<uint64_t> root_pows(N);
  root_pows[0] = 1;
  for (int i = 1; i < N; i++) {
    root_pows[i] = (root_pows[i - 1] * tables->get_root()) %
                   params.coeff_modulus()[0].value();
  }
  R rs(polytools::SealPoly(context, root_pows, &(context.first_parms_id())));
  R row(rs);
  ringsnark::linear_combination<R> sum = vars[0];
  for (int i = 1; i < N; i++) {
    sum = sum + row * vars[i];
    row *= rs;
  }
  pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(sum, 1, vars[N]));

  // Set values
  {
    BatchEncoder encoder(context);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, secret_key);
    Evaluator evaluator(context);

    vector<uint64_t> vs(N);
    Plaintext ptxt;
    Ciphertext ctxt0, ctxt1, ctxt2;

    vs[0] = 6;
    encoder.encode(vs, ptxt);
    polytools::SealPoly poly(context, ptxt, &(context.first_parms_id()));
    for (int i = 0; i < N; i++) {
      pb.val(vars[i]) = R(poly.get_coefficient_rns(i)[0]);
    }

    poly.ntt_inplace(tables);
    pb.val(vars[N]) = R(poly);
  }
  cout << "#Inputs\t" << pb.num_inputs() << endl;
  cout << "#Variables\t" << pb.num_variables() << endl;
  cout << "#Constraints\t" << pb.num_constraints() << endl;
  cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl;
  cout << endl;

  return pb;
}

ringsnark::protoboard<R> *protoboard;
ringsnark::rinocchio::keypair<R, E> *keypair;
ringsnark::rinocchio::proof<R, E> *proof;

static void BM_PlaintextCheck_Setup(benchmark::State &state) {
  ringsnark::protoboard<R> pb;
  if (!protoboard) {
    protoboard = new ringsnark::protoboard<R>(init(4096));
  }
  pb = *protoboard;

  for (auto _ : state) {
    auto kp = ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
    if (!keypair) {
      keypair = new ringsnark::rinocchio::keypair<R, E>(kp);
    }
  }
}

static void BM_PlaintextCheck_Prove(benchmark::State &state) {
  const auto pb = *protoboard;
  const auto kp = *keypair;
  for (auto _ : state) {
    auto pi = ringsnark::rinocchio::prover(kp.pk, pb.primary_input(),
                                           pb.auxiliary_input());
    if (!proof) {
      proof = new ringsnark::rinocchio::proof<R, E>(pi);
    }
  }
}

static void BM_PlaintextCheck_Verify(benchmark::State &state) {
  const auto pb = *protoboard;
  const auto kp = *keypair;
  const auto pi = *proof;
  for (auto _ : state) {
    const bool verif =
        ringsnark::rinocchio::verifier(kp.vk, pb.primary_input(), pi);
  }
}

BENCHMARK(BM_PlaintextCheck_Setup);
BENCHMARK(BM_PlaintextCheck_Prove);
BENCHMARK(BM_PlaintextCheck_Verify);
BENCHMARK_MAIN();