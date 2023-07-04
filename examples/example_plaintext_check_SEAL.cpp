#include <iostream>
#include <ringsnark/gadgetlib/protoboard.hpp>
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/seal/seal_util.hpp>
#include <ringsnark/zk_proof_systems/groth16/groth16.hpp>
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

int main() {
  EncryptionParameters params(scheme_type::bgv);
  const size_t poly_modulus_degree = 2048;
  const size_t inner_poly_modulus_degree = 8 * poly_modulus_degree;
  const size_t logT = 16;

  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_coeff_modulus(default_double_batching_modulus(
      poly_modulus_degree, inner_poly_modulus_degree));
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, logT));

  print_params(params);

  SEALContext context(params);

  typedef ringsnark::seal::RingElem R;
  typedef ringsnark::seal::EncodingElem E;

#define USE_MODSWITCH_IN_INNER_PRODUCT
  R::set_context(context);
  E::set_context(inner_poly_modulus_degree);

  const size_t num_constraints = logT + 1;
  const size_t num_variables = logT + 1;
  ringsnark::protoboard<R> pb;

  // Set public values
  ringsnark::pb_variable_array<R> vars(num_variables,
                                       ringsnark::pb_variable<R>());
  vars.allocate(pb, num_variables, "x");
  pb.set_input_sizes(1);  // plaintext is public, bit-decomposition is private

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

    vector<uint64_t> vs(poly_modulus_degree);
    Plaintext ptxt;

    vs[0] = 2;
    encoder.encode(vs, ptxt);
    auto poly = polytools::SealPoly(context, ptxt, &(context.first_parms_id()));
    vector<uint64_t> x_coeffs(ptxt.data(), ptxt.data() + poly_modulus_degree);
    // This convoluted way of initializing x is needed to trick polytools into
    // thinking that the polynomial is in NTT form (it's not).
    polytools::SealPoly x(context, x_coeffs, &(context.first_parms_id()));
    pb.val(vars[num_variables - 1]) = R(x);

    vector<uint64_t> bits(poly_modulus_degree);
    for (size_t i = 0; i < logT; i++) {
      for (size_t j = 0; j < poly_modulus_degree; j++) {
        bits[j] = (x_coeffs[j] >> i) & 1;
      }
      const polytools::SealPoly bits_poly(context, bits,
                                          &(context.first_parms_id()));
      pb.val(vars[i]) = R(bits_poly);
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