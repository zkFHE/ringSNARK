#ifndef RINGSNARK_SEAL_UTIL_HPP
#define RINGSNARK_SEAL_UTIL_HPP

#include <iostream>
#include "seal/seal.h"

using namespace std;
using namespace seal;

void print_params(const EncryptionParameters parms) {
  const SEALContext context(parms);
  cout << "[PARAM] Parameter validation (" <<  context.parameter_error_name() << "): "
       << context.parameter_error_message() << endl;
  auto qualifiers = context.first_context_data()->qualifiers();
  cout << "[PARAM] Batching enabled: " << boolalpha << qualifiers.using_batching
       << endl;
  cout << "[PARAM] poly_modulus_degree N=" << parms.poly_modulus_degree() << endl;
  cout << "[PARAM] plain_modulus log(t)=" << parms.plain_modulus().bit_count() << " bits" << endl;
  size_t coeff_modulus_bit_count = 0;
  cout << "[PARAM] coeff_modulus log(q)=";
  for (auto q_i: parms.coeff_modulus()) {
    coeff_modulus_bit_count += q_i.bit_count();
    if (q_i != parms.coeff_modulus()[0]) {
      cout << " + ";
    }
    cout << q_i.bit_count();
  }
  cout << " = " << coeff_modulus_bit_count << " bits" << endl;
  cout << "[PARAM] coeff_modulus.size=" << parms.coeff_modulus().size() << endl;

  cout << "[PARAM] plain_modulus t=" << parms.plain_modulus().value() << endl;
  cout << "[PARAM] coeff_modulus q=";
  for (auto q_i: parms.coeff_modulus()) {
    if (q_i != parms.coeff_modulus()[0]) {
      cout << " * ";
    }
    cout << q_i.value();
  }
  cout << endl;

}

#endif // RINGSNARK_SEAL_UTIL_HPP
