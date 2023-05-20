#include <string>
#include <set>

namespace ringsnark::seal_int {
[[nodiscard]] size_t RingElem::size_in_bits() const {
  return 8 * sizeof(uint64_t) * value.size();
}

bool RingElem::is_zero() const {
  for (const auto &x : value) {
    if (x != 0) {
      return false;
    }
  }
  return true;
}

bool RingElem::fast_is_zero() const { return is_zero(); }

void RingElem::negate_inplace() {
  for (size_t i = 0; i < modulus.size(); i++) {
    value[i] = seal::util::negate_uint_mod(value[i], modulus[i]);
  }
}

bool RingElem::is_invertible() const noexcept {
  uint64_t x;
  for (size_t i = 0; i < modulus.size(); i++) {
    if (!seal::util::try_invert_uint_mod(value[i], modulus[i], x)) {
      return false;
    }
  }
  return true;
}

void RingElem::invert_inplace() {
  uint64_t x;
  for (size_t i = 0; i < modulus.size(); i++) {
    if (!seal::util::try_invert_uint_mod(value[i], modulus[i], x)) {
      throw std::invalid_argument("element is not invertible in ring");
    }
    value[i] = x;
  }
}

RingElem &RingElem::operator+=(const RingElem &other) {
  for (size_t i = 0; i < modulus.size(); i++) {
    value[i] = seal::util::add_uint_mod(value[i], other.value[i], modulus[i]);
  }
  return *this;
}

RingElem &RingElem::operator-=(const RingElem &other) {
  for (size_t i = 0; i < modulus.size(); i++) {
    value[i] = seal::util::sub_uint_mod(value[i], other.value[i], modulus[i]);
  }
  return *this;
}

RingElem &RingElem::operator*=(const RingElem &other) {
  for (size_t i = 0; i < modulus.size(); i++) {
    value[i] =
        seal::util::multiply_uint_mod(value[i], other.value[i], modulus[i]);
  }
  return *this;
}

bool operator==(const RingElem &lhs, const RingElem &rhs) { return lhs.value == rhs.value; }

size_t RingElem::hash() const { throw invalid_ring_elem_types(); }

std::ostream &operator<<(std::ostream &out, const RingElem &elem) {
  out << elem.value << endl;
}

size_t EncodingElem::size_in_bits() const {
  size_t size = 0;
  for (size_t i = 0; i < ciphertexts.size(); i++) {
    auto c = ciphertexts[i];
    auto params = get_contexts()[i].get_context_data(c.parms_id())->parms();
    for (const auto &q_i : params.coeff_modulus()) {
      size += params.poly_modulus_degree() * c.size() * q_i.bit_count();
    }
  }
  return size;
}

    EncodingElem& EncodingElem::modswitch_inplace() {
        for (size_t i = 0; i < ciphertexts.size(); i++) {
            evaluators[i]->mod_switch_to_next_inplace(ciphertexts[i]);
        }
    }

std::vector<EncodingElem>
EncodingElem::encode(const SecretKey &sk, const std::vector<RingElem> &rs) {
  assert(get_contexts().size() == sk.size());
  std::vector<EncodingElem> encs(rs.size());

  // Reuse encryptors for all elements
  vector<::seal::Encryptor *> encryptors;
  encryptors.reserve(get_contexts().size());
  for (size_t i = 0; i < get_contexts().size(); i++) {
    encryptors.push_back(new ::seal::Encryptor(get_contexts()[i], sk[i]));
  }

  vector<uint64_t> vs(
      get_contexts()[0].first_context_data()->parms().poly_modulus_degree());

#pragma omp parallel for default(none)                                         \
    shared(vs, rs, encs, encryptors, ::seal::parms_id_zero)
  for (int i = 0; i < rs.size(); i++) {
    const auto &r = rs[i];
    ::seal::Plaintext ptxt;
    std::vector<::seal::Ciphertext> ctxts(get_contexts().size());

    for (size_t j = 0; j < get_contexts().size(); j++) {
      vs[0] = r.value[j];
      encoders[j]->encode(vs, ptxt);
      encryptors[j]->encrypt_symmetric(ptxt, ctxts[j]);
    }

    encs[i] = EncodingElem(ctxts);
  }
  return encs;
}

EncodingElem EncodingElem::inner_product(
    vector<EncodingElem>::const_iterator a_start,
    vector<EncodingElem>::const_iterator a_end,
    vector<RingElem>::const_iterator b_start,
    vector<RingElem>::const_iterator b_end
) {
  size_t num_moduli = EncodingElem::get_contexts()[0].first_context_data()->parms().coeff_modulus().size();
#ifdef RINGSNARK_DEBUG
  assert(a_end - a_start > 0 && "cannot compute inner product of empty vectors");
  assert(a_end - a_start == b_end - b_start && "cannot compute inner product of vectors with mismatched sizes");
  for (size_t i = 1; i < EncodingElem::get_contexts().size(); i++) {
    assert(EncodingElem::get_contexts()[i].first_context_data()->parms().coeff_modulus().size() == num_moduli);
  }
#endif

  const size_t num_modswitches = (num_moduli > 2) ? num_moduli - 2 : 0; // Number of modswitches during the sum
  set<size_t> checkpoints;
  size_t max_num_its = 0;
  vector<bool> is_zero(a_end - a_start);
  for (int i = 0; i < is_zero.size(); ++i) {
    is_zero[i] = (b_start+i)->is_zero();
    if (!is_zero[i]) {
      max_num_its++;
    }
  }

  if (num_modswitches == 1){
    checkpoints.insert(max_num_its / 2);
  } else if (num_modswitches > 1) {
    size_t delta = (a_end - a_start) / num_modswitches;
    for(int i=0; i < num_modswitches-1; ++i)
    {
      checkpoints.insert(delta * (i+1));
    }
    checkpoints.insert(max_num_its);
  }
  assert(checkpoints.size() <= num_modswitches);
  if (max_num_its > num_modswitches) {
    assert(checkpoints.size() == num_modswitches);
  }

  EncodingElem res;
  EncodingElem tmp;
  size_t curr_it_count = 0;
  for (int i = 0; i < (a_end-a_start); i++) {
    if (!is_zero[i]) {
      tmp = (*(a_start+i)) * (*(b_start+i));
#ifdef USE_MODSWITCH_IN_INNER_PRODUCT
      if (num_moduli > 1) {
        tmp.modswitch_inplace();
      }
#endif
      res += tmp;
#ifdef USE_MODSWITCH_IN_INNER_PRODUCT
      curr_it_count++;
      if (checkpoints.count(curr_it_count) == 1) {
        res.modswitch_inplace();
      }
# endif
    }
  }
  return res;
}

RingElem EncodingElem::decode(const SecretKey &sk, const EncodingElem &e) {
  // TODO: optimize to reuse same decryptor object for many invocations
  ::seal::Plaintext ptxt;
  auto parms = RingElem::get_context().first_context_data()->parms();
  std::vector<uint64_t> coeffs;

  assert(e.ciphertexts.size() == get_contexts().size());
  for (size_t i = 0; i < e.ciphertexts.size(); i++) {
    ::seal::Decryptor decryptor(get_contexts()[i], sk[i]);
    vector<uint64_t> curr_limb(parms.poly_modulus_degree());

    try {
      if (decryptor.invariant_noise_budget(e.ciphertexts[i]) <= 0) {
        // This indicates that either the parameters of the encoding scheme were
        // set to be too small, or that the prover used more budget (i.e.,
        // performed more operations) than required.
        throw decoding_error();
      }
      decryptor.decrypt(e.ciphertexts[i], ptxt);
      encoders[i]->decode(ptxt, curr_limb);
    } catch (std::invalid_argument &invalid_arg) {
      if (std::string(invalid_arg.what()) == "encrypted is empty") {
        // TODO: can we handle this case more explicitly to prevent confusion,
        // e.g., have a dedicated
        //  flag/subclass for the "0" ciphertext?
        // This should only really be an issue when the SNARK is used in non-ZK
        // mode; with ZK, the noise terms prevent the ctxt from being zero
        // w.h.p. Do nothing, curr_limb already holds all zeros.
      } else {
        throw invalid_arg;
      }
    }

    coeffs.push_back(curr_limb[0]);
  }
  assert(coeffs.size() == parms.coeff_modulus().size());
  return RingElem(coeffs);
}

EncodingElem &EncodingElem::operator+=(const EncodingElem &other) {
  if (other.is_empty()) {
    return *this;
  }
  if (this->is_empty()) {
    this->ciphertexts = other.ciphertexts;
    return *this;
  }
  assert(this->ciphertexts.size() == other.ciphertexts.size());
  assert(this->ciphertexts.size() == get_contexts().size());
  for (size_t i = 0; i < this->ciphertexts.size(); i++) {
    const bool was_ntt = this->ciphertexts[i].is_ntt_form();
    try {
      evaluators[i]->add_inplace(this->ciphertexts[i], other.ciphertexts[i]);
    } catch (std::logic_error &e) {
      if (std::string(e.what()) == "result ciphertext is transparent") {
        // Explicitly assign a zero-ciphertext
        this->ciphertexts[i] = ::seal::Ciphertext(
            get_contexts()[i], get_contexts()[i].first_parms_id());
        this->ciphertexts[i].is_ntt_form() = was_ntt;
      } else {
        throw e;
      }
    }
  }
  return *this;
}

EncodingElem &EncodingElem::operator*=(const RingElem &r) {
  ::seal::Plaintext ptxt;

  // Handle zero case explicitly, to prevent a "transparent ciphertext"
  // logic_error from SEAL
  if (r.is_zero()) {
    ciphertexts.resize(get_contexts().size());
    for (size_t i = 0; i < get_contexts().size(); i++) {
      auto context = get_contexts()[i];
      const bool was_ntt = ciphertexts[i].is_ntt_form();
      ciphertexts[i] = ::seal::Ciphertext(context, context.first_parms_id());
      ciphertexts[i].is_ntt_form() = was_ntt;
    }
    return *this;
  }

  if (r == 1) {
    return *this;
  }
  vector<uint64_t> vs(1);

  for (size_t i = 0; i < this->ciphertexts.size(); i++) {
    vs[0] = r.value[i];
    encoders[i]->encode(vs, ptxt);
    try {
      evaluators[i]->multiply_plain_inplace(this->ciphertexts[i], ptxt);
    } catch (std::logic_error err) {
      if (std::string(err.what()) == "result ciphertext is transparent") {
      } else {
        throw err;
      }
    }
  }
  return *this;
}
} // namespace ringsnark::seal_int