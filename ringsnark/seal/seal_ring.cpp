#include "seal_ring.hpp"

namespace ringsnark::seal {
    RingElem::RingElem() : value(0), type(SCALAR) {}

    RingElem::RingElem(const RingElem &other) { // Deep copy
        this->type = other.type;
        if (this->type == SCALAR) {
            this->poly = nullptr;
            this->value = other.value;
        } else if (this->type == POLY) {
            this->poly = new polytools::SealPoly(*other.poly);
            this->value = 0;
        } else {
            throw_invalid_types();
        }
    }

    RingElem::RingElem(uint64_t value) : value(value), type(SCALAR) {}

    RingElem::RingElem(::seal::SEALContext &context, ::seal::Plaintext &plaintext,
                       const ::seal::parms_id_type *parms_id_ptr)
            : poly(new polytools::SealPoly(context, plaintext, parms_id_ptr)), type(POLY) {}

    RingElem::RingElem(const polytools::SealPoly &poly) : poly(new polytools::SealPoly(poly)), type(POLY) {}

    [[nodiscard]] size_t RingElem::size_in_bits() const {
        if (this->type == SCALAR) {
            return 64;// TODO: or log2(value)?
        } else if (this->type == POLY) {
            size_t size = 0;
            for (const auto &q_i: this->poly->get_coeff_modulus()) {
                size += q_i.bit_count() * this->poly->get_coeff_count();
            }
            return size;
        } else {
            throw_invalid_types();
        }
    }

    bool RingElem::is_zero() const {
        if (this->type == SCALAR) {
            return this->value == 0;
        } else if (this->type == POLY) {
            return this->poly->is_zero();
        } else {
            throw_invalid_types();
        }
    }


    void RingElem::negate_inplace() {
        if (this->type == SCALAR) {
            // TODO: make more efficient
            this->to_poly_inplace();
            this->poly->negate_inplace();
        } else if (this->type == POLY) {
            this->poly->negate_inplace();
        } else {
            throw_invalid_types();
        }
    }


    void RingElem::invert_inplace() {
        if (this->type == SCALAR) {
            // TODO: make more efficient
            this->to_poly_inplace();
            bool success = this->poly->invert_inplace();
            if (!success) {
                throw std::invalid_argument("element is not invertible in ring");
            }
        } else if (this->type == POLY) {
            bool success = this->poly->invert_inplace();
            if (!success) {
                throw std::invalid_argument("element is not invertible in ring");
            }
        } else {
            throw_invalid_types();
        }
    }


    RingElem &RingElem::operator+=(const RingElem &other) {
        if (this->type == POLY) {
            if (other.type == POLY) {
                this->poly->add_inplace(*other.poly);
            } else if (other.type == SCALAR) {
                this->poly->add_scalar_inplace(other.value);
            } else {
                throw_invalid_types();
            }
        } else if (this->type == SCALAR) {
            if (other.type == POLY) {
                this->poly = new polytools::SealPoly(*other.poly);
                this->poly->add_scalar_inplace(this->value);
                this->type = POLY;
            } else if (other.type == SCALAR) {
                this->value += other.value;
            } else {
                throw_invalid_types();
            }
        } else {
            throw_invalid_types();
        }
        return *this;
    }


    RingElem &RingElem::operator-=(const RingElem &other) {
        if (this->type == POLY) {
            if (other.type == POLY) {
                this->poly->subtract_inplace(*other.poly);
            } else if (other.type == SCALAR) {
                this->poly->subtract_scalar_inplace(other.value);
            } else {
                throw_invalid_types();
            }
        } else if (this->type == SCALAR) {
            if (other.type == POLY) {
                this->poly = new polytools::SealPoly(*other.poly);
                this->poly->negate_inplace();
                this->poly->add_scalar_inplace(this->value);
                this->type = POLY;
            } else if (other.type == SCALAR) {
                // TODO: optimize
                this->to_poly_inplace();
                this->poly->subtract_scalar_inplace(other.value);
            } else {
                throw_invalid_types();
            }
        } else {
            throw_invalid_types();
        }
        return *this;
    }


    RingElem &RingElem::operator*=(const RingElem &other) {
        if (this->type == POLY) {
            if (other.type == POLY) {
                this->poly->multiply_inplace(*other.poly);
            } else if (other.type == SCALAR) {
                this->poly->multiply_scalar_inplace(other.value);
            } else {
                throw_invalid_types();
            }
        } else if (this->type == SCALAR) {
            if (other.type == POLY) {
                this->poly = new polytools::SealPoly(*other.poly);
                this->poly->multiply_scalar_inplace(this->value);
                this->type = POLY;
            } else if (other.type == SCALAR) {
                this->value *= other.value;
            } else {
                throw_invalid_types();
            }
        } else {
            throw_invalid_types();
        }
        return *this;
    }

    bool operator==(const RingElem &lhs, const RingElem &rhs) {
        if (lhs.type == rhs.type) {
            if (lhs.type == SCALAR) {
                return lhs.value == rhs.value;
            } else if (lhs.type == POLY) {
                return lhs.poly->is_equal(*rhs.poly);
            } else {
                RingElem::throw_invalid_types();
            }
        } else {
            // TODO: cast down to scalar instead
            if (lhs.type == SCALAR && rhs.type == POLY) {
                RingElem lhs_poly = lhs.to_poly();
                return lhs_poly.poly->is_equal(*rhs.poly);
            } else if (lhs.type == POLY && rhs.type == SCALAR) {
                RingElem rhs_poly = rhs.to_poly();
                return lhs.poly->is_equal(*rhs_poly.poly);
            } else {
                RingElem::throw_invalid_types();
            }
        }
    }


    RingElem &RingElem::to_poly_inplace() {
        if (this->type == SCALAR) {
            this->poly = new polytools::SealPoly(get_context());
            this->poly->add_scalar_inplace(this->value);
            this->type = POLY;
            assert(this->poly->is_ntt_form());
        } else if (this->type == POLY) {
            return *this;
        } else {
            throw_invalid_types();
        }
    }

    RingElem RingElem::to_scalar() const {
        if (this->type == SCALAR) {
            return *(new RingElem(*this));
        } else if (this->type == POLY) {
            // TODO: Surely there's a more efficient and simpler way to do this
            // Use SEAL's encoder to decode the current polynomial as a scalar value
            ::seal::BatchEncoder encoder(get_context());
            auto context_data = get_context().get_context_data(this->poly->get_parms_id());

            RingElem *res = new RingElem(*this);

            // Transform to iNTT form
            auto tables = context_data->small_ntt_tables();
            res->poly->intt_inplace(tables);

            std::vector<uint64_t> values(context_data->parms().poly_modulus_degree());
            ::seal::Plaintext ptxt = polytools::poly_to_ptxt(get_context(), *res->poly);
            encoder.decode(ptxt, values);

            for (int i = 1; i < values.size(); i++) {
                if (values[i] != 0) {
                    throw std::invalid_argument(
                            "cannot convert poly to scalar, coefficient at index " + std::to_string(i) + " is " +
                            std::to_string(values[i]) + " != 0");
                }
            }
            res->value = values[0];
            res->type = SCALAR;
            return *res;
        } else {
            throw_invalid_types();
        }
    }

    size_t RingElem::hash() const {
        if (this->type == SCALAR) {
            return SCALAR ^ this->value;
        } else if (this->type == POLY) {
            size_t hash = POLY;
            for (size_t i = 0; i < this->poly->get_coeff_count(); i++) {
                for (const auto &v: this->poly->get_coeff(i)) {
                    hash ^= v;
                }
            }
            return hash;
        } else {
            throw_invalid_types();
        }
    }

    std::ostream &operator<<(std::ostream &out, const RingElem &elem) {
        if (elem.type == SCALAR) {
            return out << elem.value;
        } else if (elem.type == POLY) {
            return out << elem.poly->to_json();
        } else {
            RingElem::throw_invalid_types();
        }
    }

    size_t EncodingElem::size_in_bits() const {
        size_t size = 0;
        for (const auto &c: ciphertexts) {
            auto params = context->get_context_data(c.parms_id())->parms();
            for (const auto &q_i: params.coeff_modulus()) {
                size += params.poly_modulus_degree() * c.size() * q_i.bit_count();
            }
        }
        return size;
    }

    std::vector<EncodingElem> EncodingElem::encode(const SecretKey &sk, const std::vector<RingElem> &rs) {
        std::vector<EncodingElem> encs;
        const ::seal::Encryptor encryptor(get_context(), sk); // Use secret-key encryptor for efficiency
        auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();

        std::vector<uint64_t> values(parms.poly_modulus_degree());

        ::seal::Plaintext ptxt;

        std::vector<::seal::Ciphertext> ciphertexts;
        ciphertexts.reserve(parms.coeff_modulus().size());

        for (const auto &r: rs) {
            if (r.type == SCALAR) {
                std::fill(values.begin(), values.end(), 0);
                values[0] = r.value;

                encoder->encode(values, ptxt);

                ciphertexts.resize(1);
                encryptor.encrypt_symmetric(ptxt, ciphertexts[0]);

                encs.emplace_back(ciphertexts);
            } else if (r.type == POLY) {
                // TODO: handle case where number of moduli differs, e.g., after mod-switching on the ring
                assert(r.poly->get_coeff_modulus_count() == parms.coeff_modulus().size());
                ciphertexts.resize(r.poly->get_coeff_modulus_count());
                for (size_t i = 0; i < r.poly->get_coeff_modulus_count(); i++) {
                    encoder->encode(r.poly->get_limb(i), ptxt);
                    encryptor.encrypt_symmetric(ptxt, ciphertexts[i]);
                }

                encs.emplace_back(ciphertexts);
            } else {
                RingElem::throw_invalid_types();
            }
        }

        return encs;
    }

    RingElem EncodingElem::decode(const SecretKey &sk, const EncodingElem &e) {
        // TODO: optimize to reuse same decryptor object for many invocations
        ::seal::Decryptor decryptor(get_context(), sk);
        ::seal::Plaintext ptxt;
        auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();

        assert(e.ciphertexts.size() == 1 || e.ciphertexts.size() == parms.coeff_modulus().size());
        if (e.ciphertexts.size() == 1) {
            decryptor.decrypt(e.ciphertexts[0], ptxt);

            std::vector<uint64_t> values(parms.poly_modulus_degree());
            encoder->decode(ptxt, values);

            for (size_t i = 1; i < values.size(); i++) {
                assert(values[i] == 0 && "cannot decode element to scalar ring element");
            }
            return RingElem(values[0]);
        } else {
            std::vector<uint64_t> coeffs(parms.poly_modulus_degree() * parms.coeff_modulus().size());
            for (size_t i = 0; i < e.ciphertexts.size(); i++) {
                decryptor.decrypt(e.ciphertexts[i], ptxt);
                encoder->decode(ptxt, gsl::span(coeffs.data() + i * parms.poly_modulus_degree(),
                                                parms.poly_modulus_degree()));
            }
            return RingElem(polytools::SealPoly(get_context(), coeffs, parms.parms_id()));
        }
    }

    EncodingElem &EncodingElem::operator+=(const EncodingElem &other) {
        // TODO: handle case where number of ciphertexts differ (and are > 1), e.g., in mod-switching cases?
        assert(this->ciphertexts.size() == other.ciphertexts.size() || this->ciphertexts.size() == 1 ||
               other.ciphertexts.size() == 1);
        for (size_t i = 0; i < this->ciphertexts.size(); i++) {
            evaluator->add_inplace(this->ciphertexts[i], other.ciphertexts[i]);
        }
        if (this->ciphertexts.size() < other.ciphertexts.size()) {
            // Copy ciphertexts from other in higher moduli
            for (size_t i = this->ciphertexts.size(); i < other.ciphertexts.size(); i++) {
                this->ciphertexts.emplace_back(other.ciphertexts[i]);
            }
        }
        return *this;
    }

    EncodingElem &EncodingElem::operator*=(const RingElem &r) {
        auto parms = get_context().get_context_data(get_context().first_parms_id())->parms();
        ::seal::Plaintext ptxt;

        if (r.type == SCALAR) {
            std::vector<uint64_t> values(parms.poly_modulus_degree());
            values[0] = r.value;

            encoder->encode(values, ptxt);
            for (auto &c: this->ciphertexts) {
                evaluator->multiply_plain_inplace(c, ptxt);
            }
        } else if (r.type == POLY) {
            assert(r.poly->get_coeff_modulus_count() == this->ciphertexts.size());
            for (size_t i = 0; i < this->ciphertexts.size(); i++) {
                encoder->encode(r.poly->get_limb(i), ptxt);
                evaluator->multiply_plain_inplace(this->ciphertexts[i], ptxt);
            }
        } else {
            RingElem::throw_invalid_types();
        }
        return *this;
    }


}