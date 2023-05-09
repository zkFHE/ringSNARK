#include <iostream>

#include "openfhe.h"
#include <ringsnark/openfhe/openfhe_ring.hpp>

#include <ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp>
#include <ringsnark/zk_proof_systems/groth16/groth16.hpp>

#include "poly_arith.h"
#include <ringsnark/seal/seal_ring.hpp>
#include <ringsnark/gadgetlib/protoboard.hpp>

using namespace std;
using namespace lbcrypto;

int main() {
  CCParams<CryptoContextBGVRNS> parameters;
  parameters.SetMultiplicativeDepth(3);
  parameters.SetPlaintextModulus(536903681);
  parameters.SetMaxRelinSkDeg(3);

  CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
  // enable features that you wish to use
  cryptoContext->Enable(PKE);
  cryptoContext->Enable(KEYSWITCH);
  cryptoContext->Enable(LEVELEDSHE);
  cryptoContext->Enable(ADVANCEDSHE);

  std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
  std::cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
            << std::endl;
  std::cout << "log2 q = "
            << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
            << std::endl;
}