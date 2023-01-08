#include <iostream>

#include "seal/seal.h"
#include "ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp"
#include "poly_arith.h"
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/gadgetlib/protoboard.hpp>

using namespace std;
using namespace seal;

int main() {
    EncryptionParameters params(scheme_type::bgv);
    auto poly_modulus_degree = (size_t) pow(2, 12);
    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(params);

    typedef ringsnark::seal::RingElem R;
    typedef ringsnark::seal::EncodingElem E;

    R::set_context(context);
    E::set_context();

    const size_t N = context.get_context_data(context.first_parms_id())->parms().poly_modulus_degree();
    const size_t n = 6;
    ringsnark::protoboard<R> pb;

    // Set public values
    ringsnark::pb_variable_array<R> vars(n, ringsnark::pb_variable<R>());
    vars.allocate(pb, n, "x");
    pb.set_input_sizes(n - 1); // vars[4] is private, all other values are public

    // Set constraints
    // Inputs:  x0, x1, x2, x3
    // Outputs: x4
    // Private: x5
    // x5 := x2 * x3
    // x4 := (x0 + x1) * x5
    pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(vars[2], vars[3], vars[5]));
    pb.add_r1cs_constraint(ringsnark::r1cs_constraint<R>(vars[0] + vars[1], vars[5], vars[4]));

    // Set values
    vector<ringsnark::seal::RingElem> values(n);
    {
        auto encoder = BatchEncoder(context);
        auto tables = context.get_context_data(context.first_parms_id())->small_ntt_tables();

        vector<uint64_t> vs(N);
        vector<Plaintext> xs(n);

        // Inputs
        vs[0] = 2;
        encoder.encode(vs, xs[0]);
        auto poly = polytools::SealPoly(context, xs[0], &(context.first_parms_id()));
        poly.ntt_inplace(tables);
        values[0] = ringsnark::seal::RingElem(poly);
        values[0] = ringsnark::seal::RingElem(vs[0]);

        vs[1] = 3;
        encoder.encode(vs, xs[1]);
//        poly = polytools::SealPoly(context, xs[1], &(context.first_parms_id()));
//        poly.ntt_inplace(tables);
//        values[1] = ringsnark::seal::RingElem(poly);
        values[1] = ringsnark::seal::RingElem(vs[1]);

        vs[2] = 4;
        encoder.encode(vs, xs[2]);
//        poly = polytools::SealPoly(context, xs[2], &(context.first_parms_id()));
//        poly.ntt_inplace(tables);
//        values[2] = ringsnark::seal::RingElem(poly);
        values[2] = ringsnark::seal::RingElem(vs[2]);

        vs[3] = 5;
        encoder.encode(vs, xs[3]);
//        poly = polytools::SealPoly(context, xs[3], &(context.first_parms_id()));
//        poly.ntt_inplace(tables);
//        values[3] = ringsnark::seal::RingElem(poly);
        values[3] = ringsnark::seal::RingElem(vs[3]);

        // Intermediate values
        vs[5] = vs[2] * vs[3];
//        poly = polytools::SealPoly(context, xs[2], &(context.first_parms_id()));
//        poly.multiply_inplace(polytools::SealPoly(context, xs[3], &(context.first_parms_id())));
//        values[5] = ringsnark::seal::RingElem(poly);
        values[5] = ringsnark::seal::RingElem(vs[5]);

        // Outputs
        vs[4] = (vs[0] + vs[1]) * vs[5];
//        poly = polytools::SealPoly(polytools::SealPoly(context, xs[0], &(context.first_parms_id())));
//        poly.add_inplace(polytools::SealPoly(context, xs[1], &(context.first_parms_id())));
//        poly.multiply_inplace(polytools::SealPoly(context, xs[5], &(context.first_parms_id())));
//        values[4] = ringsnark::seal::RingElem(poly);
        values[4] = ringsnark::seal::RingElem(vs[4]);
    }
    for (size_t i = 0; i < n; i++) {
        pb.val(vars[i]) = values[i];
    }
    cout << "#Inputs\t" << pb.num_inputs() << endl;
    cout << "#Variables\t" << pb.num_variables() << endl;
    cout << "#Constraints\t" << pb.num_constraints() << endl;
    assert(pb.is_satisfied());
    cout << "R1CS satisfied: " << std::boolalpha << true << endl;
    cout << endl;

    // Generate CRS
    const auto keypair = ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
    cout << "Size of pk:\t" << keypair.pk.size_in_bits() << " bits" << endl;
    cout << "Size of vk:\t" << keypair.vk.size_in_bits() << " bits" << endl;

    const auto proof = ringsnark::rinocchio::prover(keypair.pk, pb.primary_input(), pb.auxiliary_input());
    cout << "Size of proof:\t" << proof.size_in_bits() << " bits" << endl;

    const bool verif = ringsnark::rinocchio::verifier(keypair.vk, pb.primary_input(), proof);
    cout << "Verification passed: " << std::boolalpha << verif << endl;
}