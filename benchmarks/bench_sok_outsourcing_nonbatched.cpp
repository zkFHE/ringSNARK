#include <benchmark/benchmark.h>

#include <ringsnark/gadgetlib/protoboard.hpp>
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>

#include "poly_arith.h"
#include "seal/seal.h"
#include "ringsnark/seal/seal_util.hpp"

typedef ringsnark::seal::RingElem R;
typedef ringsnark::seal::EncodingElem E;
using namespace std;
using namespace seal;
using namespace ringsnark;


#define USE_MODSWITCH_IN_INNER_PRODUCT

ringsnark::protoboard<R> init(size_t N) {
    EncryptionParameters params(scheme_type::bgv);
    auto poly_modulus_degree = N;
    auto inner_poly_modulus_degree = 8 * poly_modulus_degree;

    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(default_double_batching_modulus(
            poly_modulus_degree, inner_poly_modulus_degree));
    params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 16));
    SEALContext context(params);

    print_params(params);

    try {
        R::set_context(context);
        E::set_context(inner_poly_modulus_degree);
    } catch (std::invalid_argument &e) {
    }
    const size_t num_features = 64;

    const size_t num_input_variables = 2 * num_features + 5;
    ringsnark::protoboard<R> pb;

    vector<pb_variable_array<R>> in1(num_features), in2(num_features);
    for (size_t i = 0; i < num_features; i++) {
        in1[i].allocate(pb, 2, "in" + to_string(i));
        in2[i].allocate(pb, 2, "in" + to_string(i));
    }
    pb_variable_array<R> out;
    out.allocate(pb, 5, "out");
    pb.set_input_sizes(num_input_variables);

    vector<linear_combination<R>> sum(3);
    // multiply and aggregate
    for (size_t i = 0; i < num_features; i++) {
        pb_variable<R> prod0, prod2;
        prod0.allocate(pb, "in" + to_string(i));
        prod2.allocate(pb, "in" + to_string(i));

        pb_variable_array<R> tmp;
        tmp.allocate(pb, 2, "tmp" + to_string(i));

        pb.add_r1cs_constraint(r1cs_constraint<R>(in1[i][0], in2[i][0], prod0));
        sum[0] = sum[0] + prod0;

        pb.add_r1cs_constraint(r1cs_constraint<R>(in1[i][0], in2[i][1], tmp[0]));
        pb.add_r1cs_constraint(r1cs_constraint<R>(in1[i][1], in2[i][0], tmp[1]));
        sum[1] = sum[1] + tmp[0] + tmp[1];

        pb.add_r1cs_constraint(r1cs_constraint<R>(in1[i][1], in2[i][1], prod2));
        sum[2] = sum[2] + prod2;
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

    cout << "#inputs\t" << pb.num_inputs() << endl;
    cout << "#variables\t" << pb.num_variables() << endl;
    cout << "#constraints\t" << pb.num_constraints() << endl;
    cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl;
    cout << endl;

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

//        pb.val(x0) = R(polytools::SealPoly(context, ctxt0, 0));
//        pb.val(x1) = R(polytools::SealPoly(context, ctxt0, 1));
//        pb.val(y0) = R(polytools::SealPoly(context, ctxt1, 0));
//        pb.val(y1) = R(polytools::SealPoly(context, ctxt1, 1));
//        pb.val(z0) = R(polytools::SealPoly(context, ctxt2, 0));
//        pb.val(z1) = R(polytools::SealPoly(context, ctxt2, 1));
//        pb.val(z2) = R(polytools::SealPoly(context, ctxt2, 2));
//
//        auto poly = polytools::SealPoly(context, ctxt0, 1);
//        poly.multiply_inplace(polytools::SealPoly(context, ctxt1, 0));
//        pb.val(tmp) = R(poly);
    }
    cout << "#Inputs\t" << pb.num_inputs() << endl;
    cout << "#Variables\t" << pb.num_variables() << endl;
    cout << "#Constraints\t" << pb.num_constraints() << endl;
    cout << "#Constraints\t" << pb.num_constraints() << endl;
    cout << "R1CS satisfied: " << std::boolalpha << pb.is_satisfied() << endl;
    cout << endl;

    return pb;
}

ringsnark::protoboard<R> *global_protoboard;
ringsnark::rinocchio::keypair<R, E> *global_keypair;
ringsnark::rinocchio::proof<R, E> *global_proof;

static void BM_PlaintextCheck_Setup(benchmark::State &state) {
    ringsnark::protoboard<R> pb;
    if (!global_protoboard) {
        global_protoboard = new ringsnark::protoboard<R>(init(2048));
    }
    pb = *global_protoboard;

    for (auto _: state) {
        auto kp = ringsnark::rinocchio::generator<R, E>(pb.get_constraint_system());
        if (!global_keypair) {
            global_keypair = new ringsnark::rinocchio::keypair<R, E>(kp);
        }
    }
}

static void BM_PlaintextCheck_Prove(benchmark::State &state) {
    const auto pb = *global_protoboard;
    const auto kp = *global_keypair;
    for (auto _: state) {
        auto pi = ringsnark::rinocchio::prover(kp.pk, pb.primary_input(),
                                               pb.auxiliary_input());
        if (!global_proof) {
            global_proof = new ringsnark::rinocchio::proof<R, E>(pi);
        }
    }
}

static void BM_PlaintextCheck_Verify(benchmark::State &state) {
    const auto pb = *global_protoboard;
    const auto kp = *global_keypair;
    const auto pi = *global_proof;
    for (auto _: state) {
        const bool verif =
                ringsnark::rinocchio::verifier(kp.vk, pb.primary_input(), pi);
    }
}

BENCHMARK(BM_PlaintextCheck_Setup);
BENCHMARK(BM_PlaintextCheck_Prove);
BENCHMARK(BM_PlaintextCheck_Verify);
BENCHMARK_MAIN();