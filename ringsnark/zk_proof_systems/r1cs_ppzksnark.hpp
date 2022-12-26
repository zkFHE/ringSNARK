/** @file
 *****************************************************************************

 Declaration of interfaces for a pre-processing zero-knowledge SNARK (ppzkSNARK) for R1CS.

 This includes:
 - class for proving key
 - class for verification key
 - class for processed verification key
 - class for key pair (proving key & verification key)
 - class for proof
 - generator algorithm
 - prover algorithm
 - verifier algorithm


 Acronyms:

 - R1CS = "Rank-1 Constraint Systems"
 - ppzkSNARK = "PreProcessing Zero-Knowledge Succinct Non-interactive ARgument of Knowledge"

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef R1CS_PPZKSNARK_HPP_
#define R1CS_PPZKSNARK_HPP_

#include <memory>

#include <ringsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>

namespace ringsnark {

/******************************** Proving key ********************************/

    template<typename RingT, typename EncT>
    class proving_key;

    template<typename RingT, typename EncT>
    std::ostream &operator<<(std::ostream &out, const proving_key<RingT, EncT> &pk);

    template<typename RingT, typename EncT>
    std::istream &operator>>(std::istream &in, proving_key<RingT, EncT> &pk);

/**
 * A proving key for the R1CS ppzkSNARK.
 */
    template<typename RingT, typename EncT>
    class proving_key {
    public:
        proving_key() = default;

        [[nodiscard]] virtual size_t size_in_bits() const = 0;
    };

    template<typename RingT, typename EncT>
    bool operator==(const proving_key<RingT, EncT> &lhs, const proving_key<RingT, EncT> &rhs);


/******************************* Verification key ****************************/

    template<typename RingT, typename EncT>
    class verification_key;

    template<typename RingT, typename EncT>
    std::ostream &operator<<(std::ostream &out, const verification_key<RingT, EncT> &vk);

    template<typename RingT, typename EncT>
    std::istream &operator>>(std::istream &in, verification_key<RingT, EncT> &vk);

/**
 * A verification key for the R1CS ppzkSNARK.
 */
    template<typename RingT, typename EncT>
    class verification_key {
    public:
        [[nodiscard]] virtual size_t size_in_bits() const = 0;
    };


/************************ Processed verification key *************************/

    template<typename RingT, typename EncT>
    class processed_verification_key;

    template<typename RingT, typename EncT>
    std::ostream &operator<<(std::ostream &out, const processed_verification_key<RingT, EncT> &pvk);

    template<typename RingT, typename EncT>
    std::istream &operator>>(std::istream &in, processed_verification_key<RingT, EncT> &pvk);

/**
 * A processed verification key for the R1CS ppzkSNARK.
 *
 * Compared to a (non-processed) verification key, a processed verification key
 * contains a small constant amount of additional pre-computed information that
 * enables a faster verification time.
 */
    template<typename RingT, typename EncT>
    class processed_verification_key {
    public:
    };


/********************************** Key pair *********************************/

/**
 * A key pair for the R1CS ppzkSNARK, which consists of a proving key and a verification key.
 */
//    template<typename RingT, typename EncT>
//    using keypair = std::tuple<proving_key<RingT, EncT>, verification_key<RingT, EncT>>;

    template<typename RingT, typename EncT>
    class keypair {
    public:
        proving_key<RingT, EncT> &pk;
        verification_key<RingT, EncT> &vk;

        keypair(const keypair<RingT, EncT> &other) = default;

        keypair(keypair<RingT, EncT> &&other) noexcept = default;

        keypair(proving_key<RingT, EncT> &pk, verification_key<RingT, EncT> &vk) : pk(pk), vk(vk) {}

        keypair(proving_key<RingT, EncT> &&pk, verification_key<RingT, EncT> &&vk) : pk(std::move(pk)),
                                                                                     vk(std::move(vk)) {}
    };


/*********************************** Proof ***********************************/

    template<typename RingT, typename EncT>
    class proof;

    template<typename RingT, typename EncT>
    std::ostream &operator<<(std::ostream &out, const proof<RingT, EncT> &proof);

    template<typename RingT, typename EncT>
    std::istream &operator>>(std::istream &in, proof <RingT, EncT> &proof);

/**
 * A proof for the R1CS ppzkSNARK.
 *
 * While the proof has a structure, externally one merely opaquely produces,
 * serializes/deserializes, and verifies proofs. We only expose some information
 * about the structure for statistics purposes.
 */
    template<typename RingT, typename EncT>
    class proof {
    public:

        [[nodiscard]] virtual size_t size_in_bits() const = 0;

//        [[nodiscard]] bool is_well_formed() const = 0;

//        bool operator==(const proof<RingT, EncT> &other) const;
    };


/***************************** Main algorithms *******************************/

/**
 * A generator algorithm for the R1CS ppzkSNARK.
 *
 * Given a R1CS constraint system CS, this algorithm produces proving and verification keys for CS.
 */
    template<typename RingT, typename EncT>
    keypair<RingT, EncT> generator(const r1cs_constraint_system<RingT> &cs);

/**
 * A prover algorithm for the R1CS ppzkSNARK.
 *
 * Given a R1CS primary input X and a R1CS auxiliary input Y, this algorithm
 * produces a proof (of knowledge) that attests to the following statement:
 *               ``there exists Y such that CS(X,Y)=0''.
 * Above, CS is the R1CS constraint system that was given as input to the generator algorithm.
 */
    template<typename RingT, typename EncT>
    proof<RingT, EncT> prover(const proving_key<RingT, EncT> &pk,
                              const r1cs_primary_input<RingT> &primary_input,
                              const r1cs_auxiliary_input<RingT> &auxiliary_input);

/*
 Below are four variants of verifier algorithm for the R1CS ppzkSNARK.

 These are the four cases that arise from the following two choices:

 (1) The verifier accepts a (non-processed) verification key or, instead, a processed verification key.
     In the latter case, we call the algorithm an "online verifier".

 (2) The verifier checks for "weak" input consistency or, instead, "strong" input consistency.
     Strong input consistency requires that |primary_input| = CS.num_inputs, whereas
     weak input consistency requires that |primary_input| <= CS.num_inputs (and
     the primary input is implicitly padded with zeros up to length CS.num_inputs).
 */

/**
 * A verifier algorithm for the R1CS ppzkSNARK that:
 * (1) accepts a non-processed verification key, and
 * (2) has weak input consistency.
 */
    template<typename RingT, typename EncT>
    bool verifier_weak_IC(const verification_key<RingT, EncT> &vk,
                          const r1cs_primary_input<RingT> &primary_input,
                          const proof<RingT, EncT> &proof);

/**
 * A verifier algorithm for the R1CS ppzkSNARK that:
 * (1) accepts a non-processed verification key, and
 * (2) has strong input consistency.
 */
    template<typename RingT, typename EncT>
    bool verifier_strong_IC(const verification_key<RingT, EncT> &vk,
                            const r1cs_primary_input<RingT> &primary_input,
                            const proof<RingT, EncT> &proof);

/**
 * Convert a (non-processed) verification key into a processed verification key.
 */
    template<typename RingT, typename EncT>
    processed_verification_key<RingT, EncT> verifier_process_vk(const verification_key<RingT, EncT> &vk);

/**
 * A verifier algorithm for the R1CS ppzkSNARK that:
 * (1) accepts a processed verification key, and
 * (2) has weak input consistency.
 */
    template<typename RingT, typename EncT>
    bool online_verifier_weak_IC(const processed_verification_key<RingT, EncT> &pvk,
                                 const r1cs_primary_input<RingT> &input,
                                 const proof<RingT, EncT> &proof);

/**
 * A verifier algorithm for the R1CS ppzkSNARK that:
 * (1) accepts a processed verification key, and
 * (2) has strong input consistency.
 */
    template<typename RingT, typename EncT>
    bool online_verifier_strong_IC(const processed_verification_key<RingT, EncT> &pvk,
                                   const r1cs_primary_input<RingT> &primary_input,
                                   const proof<RingT, EncT> &proof);


} // ringsnark

#endif // R1CS_PPZKSNARK_HPP_
