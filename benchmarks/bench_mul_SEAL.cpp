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
  params.set_plain_modulus(PlainModulus::Batching(N, 32));
  SEALContext context(params);

  try {
    R::set_context(context);
    E::set_context();
  } catch (std::invalid_argument &e) {
  }

  const size_t num_input_variables = 2 + 2 + 3;
  const size_t num_variables = 2 + 2 + 3 + 1;
  ringsnark::protoboard<R> pb;

  // Set public values
  ringsnark::pb_variable_array<R> vars(num_variables,
                                       ringsnark::pb_variable<R>());
  vars.allocate(pb, num_variables, "x");
  pb.set_input_sizes(num_input_variables);
  // For (z0, z1, z2) = (x0, x1) * (y0, y1), the variables are arranged as
  // follows: [x0, x1, y0, y1, tmp, z0, z1, z2], with z0 = x0 * y0 tmp = x0 * y1
  // z1 = tmp + (x1 * y0)
  // z2 = x1 * y1
  auto x0 = vars[0], x1 = vars[1], y0 = vars[2], y1 = vars[3], tmp = vars[4],
       z0 = vars[5], z1 = vars[6], z2 = vars[7];

  // Set constraints
  pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(x0, y0, z0));
  pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(x1, y0, tmp));
  pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(x0, y1, z1 - tmp));
  pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(x1, y1, z2));

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
    encryptor.encrypt_symmetric(ptxt, ctxt0);
    assert(ctxt0.is_ntt_form());

    vs[0] = 7;
    encoder.encode(vs, ptxt);
    encryptor.encrypt_symmetric(ptxt, ctxt1);
    assert(ctxt1.is_ntt_form());

    evaluator.multiply(ctxt0, ctxt1, ctxt2);

    pb.val(x0) = R(polytools::SealPoly(context, ctxt0, 0));
    pb.val(x1) = R(polytools::SealPoly(context, ctxt0, 1));
    pb.val(y0) = R(polytools::SealPoly(context, ctxt1, 0));
    pb.val(y1) = R(polytools::SealPoly(context, ctxt1, 1));
    pb.val(z0) = R(polytools::SealPoly(context, ctxt2, 0));
    pb.val(z1) = R(polytools::SealPoly(context, ctxt2, 1));
    pb.val(z2) = R(polytools::SealPoly(context, ctxt2, 2));

    auto poly = polytools::SealPoly(context, ctxt0, 1);
    poly.multiply_inplace(polytools::SealPoly(context, ctxt1, 0));
    pb.val(tmp) = R(poly);
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
    protoboard = new ringsnark::protoboard<R>(init(8192));
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