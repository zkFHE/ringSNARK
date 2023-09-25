#include <iostream>
#include <ringsnark/gadgetlib/protoboard.hpp>
#include <ringsnark/seal/seal_util.hpp>
#include <ringsnark/seal_int/seal_ring.hpp>
#include <ringsnark/zk_proof_systems/groth16/groth16.hpp>
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

int main() {
  EncryptionParameters params(scheme_type::bgv);
  auto poly_modulus_degree = 1024;
  auto inner_poly_modulus_degree = 8 * poly_modulus_degree;

  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_coeff_modulus(default_double_batching_modulus(
      poly_modulus_degree, inner_poly_modulus_degree));
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 16));

  print_params(params);

  SEALContext context(params);

  typedef ringsnark::seal_int::RingElem R;
  typedef ringsnark::seal_int::EncodingElem E;

#define USE_MODSWITCH_IN_INNER_PRODUCT
  R::set_context(context);
  E::set_context(inner_poly_modulus_degree);

  const size_t N = context.get_context_data(context.first_parms_id())
                       ->parms()
                       .poly_modulus_degree();
  const size_t n = 6;
  ringsnark::protoboard<R> pb;

  // Set public values
  ringsnark::pb_variable_array<R> vars(n, ringsnark::pb_variable<R>());
  vars.allocate(pb, n * N, "x");
  pb.set_input_sizes((n - 1) *
                     N);  // vars[4] is private, all other values are public

  // Set constraints
  // Inputs:  x0, x1, x2, x3
  // Outputs: x4
  // Private: x5
  // x5 := x2 * x3
  // x4 := (x0 + x1) * x5
  for (int i = 0; i < N; i++) {
    pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(
        vars[2 * N + i], vars[3 * N + i], vars[5 * N + i]));
    pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(
        vars[0 * N + i] + vars[1 * N + i], vars[5 * N + i], vars[4 * N + i]));
  }

  // Set values
  vector<ringsnark::seal_int::RingElem> values(n);
  {
    auto encoder = BatchEncoder(context);
    auto tables =
        context.get_context_data(context.first_parms_id())->small_ntt_tables();

    vector<uint64_t> vs(N);
    Plaintext ptxt;
    vector<::polytools::SealPoly> polys(n, ::polytools::SealPoly(context));

    // Inputs
    vs[0] = 2;
    encoder.encode(vs, ptxt);
    auto poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    poly.ntt_inplace(tables);
    polys[0] = poly;

    vs[1] = 3;
    encoder.encode(vs, ptxt);
    poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    poly.ntt_inplace(tables);
    polys[1] = poly;

    vs[2] = 4;
    encoder.encode(vs, ptxt);
    poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    poly.ntt_inplace(tables);
    polys[2] = poly;

    vs[3] = 5;
    encoder.encode(vs, ptxt);
    poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    poly.ntt_inplace(tables);
    polys[3] = poly;

    // Intermediate values
    vs[5] = vs[2] * vs[3];
    poly = ::polytools::SealPoly(polys[2]);
    poly.multiply_inplace(polys[3]);
    polys[5] = poly;

    // Outputs
    vs[4] = (vs[0] + vs[1]) * vs[5];
    poly = ::polytools::SealPoly(polys[0]);
    poly.add_inplace(polys[1]);
    poly.multiply_inplace(polys[5]);
    polys[4] = poly;

    for (size_t j = 0; j < polys.size(); j++) {
      for (size_t i = 0; i < N; i++) {
        pb.val(vars[j * N + i]) = R(polys[j].get_coefficient_rns(i));
      }
    }
  }

  cout << "#Inputs\t" << pb.num_inputs() << endl;
  cout << "#Variables\t" << pb.num_variables() << endl;
  cout << "#Constraints\t" << pb.num_constraints() << endl;
  cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl;
  cout << endl;

  {
    cout << "=== Rinocchio ===" << endl;
    const auto keypair =
        ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
    cout << "Size of pk:\t" << keypair.pk.size_in_bits() << " bits" << endl;
    cout << "Size of vk:\t" << keypair.vk.size_in_bits() << " bits" << endl;

    const auto proof = ringsnark::rinocchio::prover(
        keypair.pk, pb.primary_input(), pb.auxiliary_input());
    cout << "Size of proof:\t" << proof.size_in_bits() << " bits" << endl;

    const bool verif =
        ringsnark::rinocchio::verifier(keypair.vk, pb.primary_input(), proof);
    cout << "Verification passed: " << std::boolalpha << verif << endl;
  }
  {
    cout << "=============" << endl;
    cout << "=== RingGroth16 ===" << endl;
    const auto keypair =
        ringsnark::groth16::generator<R, E>(pb.get_constraint_system());
    cout << "Size of pk:\t" << keypair.pk.size_in_bits() << " bits" << endl;
    cout << "Size of vk:\t" << keypair.vk.size_in_bits() << " bits" << endl;

    const auto proof = ringsnark::groth16::prover(
        keypair.pk, pb.primary_input(), pb.auxiliary_input());
    cout << "Size of proof:\t" << proof.size_in_bits() << " bits" << endl;

    const bool verif =
        ringsnark::groth16::verifier(keypair.vk, pb.primary_input(), proof);
    cout << "Verification passed: " << std::boolalpha << verif << endl;
  }
}