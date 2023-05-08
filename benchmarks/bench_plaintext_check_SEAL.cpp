#include "seal/seal.h"
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/gadgetlib/protoboard.hpp>
#include <benchmark/benchmark.h>

using namespace std;
using namespace seal;

typedef ringsnark::seal::RingElem R;
typedef ringsnark::seal::EncodingElem E;

ringsnark::protoboard<R> init(size_t logT) {
  EncryptionParameters params(scheme_type::bgv);
  auto poly_modulus_degree = (size_t)pow(2, 11);
  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, logT));
  SEALContext context(params);

  try {
    R::set_context(context);
    E::set_context();
  } catch (std::invalid_argument& e) {}

//  const size_t logT = params.plain_modulus().bit_count();
  const size_t N = context.get_context_data(context.first_parms_id())
                       ->parms()
                       .poly_modulus_degree();
  const size_t num_constraints = logT + 1;
  const size_t num_variables = logT + 1;
  ringsnark::protoboard<R> pb;

  // Set public values
  ringsnark::pb_variable_array<R> vars(num_variables,
                                       ringsnark::pb_variable<R>());
  vars.allocate(pb, num_variables, "x");
  // plaintext is public, bit-decomposition is private
  pb.set_input_sizes(num_variables - 1);

  // Set constraints
  ringsnark::linear_combination<R> sum = 0;
  for (size_t i = 0; i < logT; i++) {
    pb.add_r1cs_constraint(
        ringsnark::r1cs_constraint<R>(vars[i], 1 - vars[i], 0));
    sum = sum + (1 << i) * vars[i];
  }
  pb.add_r1cs_constraint(
      ringsnark::r1cs_constraint<R>(vars[num_variables - 1], 1, sum));

  // Set values
  {
    auto encoder = BatchEncoder(context);

    vector<uint64_t> vs(N);
    Plaintext ptxt;

    vs[0] = 2;
    encoder.encode(vs, ptxt);
    auto poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    vector<uint64_t> x_coeffs(ptxt.data(), ptxt.data() + N);
    // This convoluted way of initializing x is needed to trick polytools into thinking that the polynomial is in NTT form (it's not).
    polytools::SealPoly x(context, x_coeffs, &(context.first_parms_id()));
    pb.val(vars[num_variables - 1]) = R(x);

    vector<uint64_t> bits(N);
    for (size_t i = 0; i < logT; i++) {
      for (size_t j = 0; j < N; j++) {
        bits[j] = (x_coeffs[j] >> i) & 1;
      }
      const polytools::SealPoly bits_poly(context, bits,
                                          &(context.first_parms_id()));
      pb.val(vars[i]) = R(bits_poly);
    }
  }
  //    cout << "#Inputs\t" << pb.num_inputs() << endl;
  //    cout << "#Variables\t" << pb.num_variables() << endl;
  //    cout << "#Constraints\t" << pb.num_constraints() << endl;
  //    cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl; cout << endl;
  //
  return pb;
}

ringsnark::protoboard<R>* protoboard;
ringsnark::rinocchio::keypair<R, E>* keypair;
ringsnark::rinocchio::proof<R, E>* proof;

static void BM_PlaintextCheck_Setup(benchmark::State& state) {
    ringsnark::protoboard<R> pb;
    if (!protoboard) {
        protoboard = new ringsnark::protoboard<R>(init(16));
    }
    pb = *protoboard;

    for (auto _ : state) {
        auto kp = ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
        if (!keypair) {
          keypair = new ringsnark::rinocchio::keypair<R, E>(kp);
        }
    }

}

static void BM_PlaintextCheck_Prove(benchmark::State& state) {
    const auto pb = *protoboard;
    const auto kp = *keypair;
    for (auto _ : state) {
        auto pi = ringsnark::rinocchio::prover(kp.pk, pb.primary_input(), pb.auxiliary_input());
        if (!proof) {
          proof = new ringsnark::rinocchio::proof<R, E>(pi);
        }
    }

}

static void BM_PlaintextCheck_Verify(benchmark::State& state) {
    const auto pb = *protoboard;
    const auto kp = *keypair;
    const auto pi = *proof;
    for (auto _ : state) {
        const bool verif = ringsnark::rinocchio::verifier(kp.vk, pb.primary_input(), pi);
    }
}

BENCHMARK(BM_PlaintextCheck_Setup);
BENCHMARK(BM_PlaintextCheck_Prove);
BENCHMARK(BM_PlaintextCheck_Verify);
BENCHMARK_MAIN();