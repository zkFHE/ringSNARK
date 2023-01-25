/** @file
*****************************************************************************

Implementation of interfaces for a QRP ("Quadratic Ring Program").

See qrp.hpp .

*****************************************************************************/

#ifndef QRP_TCC_
#define QRP_TCC_

namespace ringsnark {
    template<typename T, typename U>
    T inner_product(typename std::vector<T>::const_iterator a_start,
                    typename std::vector<T>::const_iterator a_end,
                    typename std::vector<U>::const_iterator b_start,
                    typename std::vector<U>::const_iterator b_end) {
        assert(a_end - a_start > 0 && "cannot compute inner product of empty vectors");
        assert(a_end - a_start == b_end - b_start && "cannot compute inner product of vectors with mismatched sizes");
        auto a_it = a_start;
        auto b_it = b_start;
        T res(*a_start);
        a_it++;
        res *= *b_start;
        b_it++;
        while (a_it != a_end) {
            res += (*a_it) * (*b_it);
            a_it++;
            b_it++;
        }
        return res;
    }

    template<typename T>
    T inner_product(typename std::vector<T>::const_iterator a_start,
                    typename std::vector<T>::const_iterator a_end,
                    typename std::vector<T>::const_iterator b_start,
                    typename std::vector<T>::const_iterator b_end) {
        return inner_product<T, T>(a_start, a_end, b_start, b_end);
    }


    template<typename RingT>
    qrp_instance<RingT>::qrp_instance(
            const std::shared_ptr<evaluation_domain<RingT>> &domain,
            const size_t num_variables,
            const size_t degree,
            const size_t num_inputs,
            const std::vector<std::map<size_t, RingT>> &A_in_Lagrange_basis,
            const std::vector<std::map<size_t, RingT>> &B_in_Lagrange_basis,
            const std::vector<std::map<size_t, RingT>> &C_in_Lagrange_basis
    ) :
            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            domain(domain),
            A_in_Lagrange_basis(A_in_Lagrange_basis),
            B_in_Lagrange_basis(B_in_Lagrange_basis),
            C_in_Lagrange_basis(C_in_Lagrange_basis) {
    }

    template<typename RingT>
    qrp_instance<RingT>::qrp_instance(const std::shared_ptr<evaluation_domain<RingT>> &domain,
                                      const size_t num_variables,
                                      const size_t degree,
                                      const size_t num_inputs,
                                      std::vector<std::map<size_t, RingT>>
                                      &&A_in_Lagrange_basis,
                                      std::vector<std::map<size_t, RingT>> &&B_in_Lagrange_basis,
                                      std::vector<std::map<size_t, RingT>>
                                      &&C_in_Lagrange_basis) :

            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            domain(domain),
            A_in_Lagrange_basis(std::move(A_in_Lagrange_basis)),
            B_in_Lagrange_basis(std::move(B_in_Lagrange_basis)),
            C_in_Lagrange_basis(std::move(C_in_Lagrange_basis)) {
    }

    template<typename RingT>
    size_t qrp_instance<RingT>::num_variables() const {
        return num_variables_;
    }

    template<typename RingT>
    size_t qrp_instance<RingT>::degree() const {
        return degree_;
    }

    template<typename RingT>
    size_t qrp_instance<RingT>::num_inputs() const {
        return num_inputs_;
    }

    template<typename RingT>
    bool qrp_instance<RingT>::is_satisfied(const qrp_witness<RingT> &witness) const {
        const RingT t = RingT::random_element();

        std::vector<RingT> At(this->num_variables() + 1, RingT::zero());
        std::vector<RingT> Bt(this->num_variables() + 1, RingT::zero());
        std::vector<RingT> Ct(this->num_variables() + 1, RingT::zero());
        std::vector<RingT> Ht(this->degree() + 1);

        const RingT Zt = this->domain->compute_vanishing_polynomial(t);

        const std::vector<RingT> u = this->domain->evaluate_all_lagrange_polynomials(t);

        for (size_t i = 0; i < this->num_variables() + 1; ++i) {
            for (auto &el: A_in_Lagrange_basis[i]) {
                At[i] += u[el.first] * el.second;
            }

            for (auto &el: B_in_Lagrange_basis[i]) {
                Bt[i] += u[el.first] * el.second;
            }

            for (auto &el: C_in_Lagrange_basis[i]) {
                Ct[i] += u[el.first] * el.second;
            }
        }

        RingT ti = RingT::one();
        for (size_t i = 0; i < this->degree() + 1; ++i) {
            Ht[i] = ti;
            ti *= t;
        }

        const qrp_instance_evaluation<RingT> eval_qrp_inst(this->domain,
                                                           this->num_variables(),
                                                           this->degree(),
                                                           this->num_inputs(),
                                                           t,
                                                           std::move(At),
                                                           std::move(Bt),
                                                           std::move(Ct),
                                                           std::move(Ht),
                                                           Zt);
        return eval_qrp_inst.is_satisfied(witness);
    }

    template<typename RingT>
    qrp_instance_evaluation<RingT>::qrp_instance_evaluation(const std::shared_ptr<evaluation_domain<RingT>> &domain,
                                                            const size_t num_variables,
                                                            const size_t degree,
                                                            const size_t num_inputs,
                                                            const RingT &t,
                                                            const std::vector<RingT> &At,
                                                            const std::vector<RingT> &Bt,
                                                            const std::vector<RingT> &Ct,
                                                            const std::vector<RingT> &Ht,
                                                            const RingT &Zt
    ) :

            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            domain(domain),
            t(t),
            At(At),
            Bt(Bt),
            Ct(Ct),
            Ht(Ht),
            Zt(Zt) {
    }

    template<typename RingT>
    qrp_instance_evaluation<RingT>::qrp_instance_evaluation(const std::shared_ptr<evaluation_domain<RingT>> &domain,
                                                            const size_t num_variables,
                                                            const size_t degree,
                                                            const size_t num_inputs,
                                                            const RingT &t,
                                                            std::vector<RingT>
                                                            &&At,
                                                            std::vector<RingT> &&Bt,
                                                            std::vector<RingT>
                                                            &&Ct,
                                                            std::vector<RingT> &&Ht,
                                                            const RingT &Zt) :
            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            domain(domain),
            t(t),
            At(std::move(At)),
            Bt(std::move(Bt)),
            Ct(std::move(Ct)),
            Ht(std::move(Ht)),
            Zt(Zt) {
    }

    template<typename RingT>
    size_t qrp_instance_evaluation<RingT>::num_variables() const {
        return num_variables_;
    }

    template<typename RingT>
    size_t qrp_instance_evaluation<RingT>::degree() const {
        return degree_;
    }

    template<typename RingT>
    size_t qrp_instance_evaluation<RingT>::num_inputs() const {
        return num_inputs_;
    }

    template<typename RingT>
    bool qrp_instance_evaluation<RingT>::is_satisfied(const qrp_witness<RingT> &witness) const {

        if (this->num_variables() != witness.num_variables()) {
            return false;
        }

        if (this->degree() != witness.degree()) {
            return false;
        }

        if (this->num_inputs() != witness.num_inputs()) {
            return false;
        }

        if (this->num_variables() != witness.coefficients_for_ABCs.size()) {
            return false;
        }

        if (this->degree() + 1 != witness.coefficients_for_H.size()) {
            return false;
        }

        if (this->At.size() != this->num_variables() + 1 || this->Bt.size() != this->num_variables() + 1 ||
            this->Ct.size() != this->num_variables() + 1) {
            return false;
        }

        if (this->Ht.size() != this->degree() + 1) {
            return false;
        }

        if (this->Zt != this->domain->compute_vanishing_polynomial(this->t)) {
            return false;
        }

        RingT ans_A = this->At[0] + witness.d1 * this->Zt;
        RingT ans_B = this->Bt[0] + witness.d2 * this->Zt;
        RingT ans_C = this->Ct[0] + witness.d3 * this->Zt;
        RingT ans_H = RingT::zero();

        ans_A = ans_A + inner_product<RingT>(this->At.begin() + 1,
                                             this->At.begin() + 1 + this->num_variables(),
                                             witness.coefficients_for_ABCs.begin(),
                                             witness.coefficients_for_ABCs.begin() + this->num_variables());
        ans_B = ans_B + inner_product<RingT>(this->Bt.begin() + 1,
                                             this->Bt.begin() + 1 + this->num_variables(),
                                             witness.coefficients_for_ABCs.begin(),
                                             witness.coefficients_for_ABCs.begin() + this->num_variables());
        ans_C = ans_C + inner_product<RingT>(this->Ct.begin() + 1,
                                             this->Ct.begin() + 1 + this->num_variables(),
                                             witness.coefficients_for_ABCs.begin(),
                                             witness.coefficients_for_ABCs.begin() + this->num_variables());
        ans_H = ans_H + inner_product<RingT>(this->Ht.begin(),
                                             this->Ht.begin() + this->degree() + 1,
                                             witness.coefficients_for_H.begin(),
                                             witness.coefficients_for_H.begin() + this->degree() + 1);

        if (ans_A * ans_B - ans_C != ans_H * this->Zt) {
            return false;
        }

        return true;
    }

    template<typename RingT>
    qrp_witness<RingT>::qrp_witness(const size_t num_variables,
                                    const size_t degree,
                                    const size_t num_inputs,
                                    const RingT &d1,
                                    const RingT &d2,
                                    const RingT &d3,
                                    const std::vector<RingT> &coefficients_for_ABCs,
                                    const std::vector<RingT> &coefficients_for_A_io,
                                    const std::vector<RingT> &coefficients_for_B_io,
                                    const std::vector<RingT> &coefficients_for_C_io,
                                    const std::vector<RingT> &coefficients_for_A_mid,
                                    const std::vector<RingT> &coefficients_for_B_mid,
                                    const std::vector<RingT> &coefficients_for_C_mid,
                                    const std::vector<RingT> &coefficients_for_Z,
                                    const std::vector<RingT> &coefficients_for_H) :
            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            d1(d1),
            d2(d2),
            d3(d3),
            coefficients_for_ABCs(coefficients_for_ABCs),
            coefficients_for_A_io(coefficients_for_A_io),
            coefficients_for_B_io(coefficients_for_B_io),
            coefficients_for_C_io(coefficients_for_C_io),
            coefficients_for_A_mid(coefficients_for_A_mid),
            coefficients_for_B_mid(coefficients_for_B_mid),
            coefficients_for_C_mid(coefficients_for_C_mid),
            coefficients_for_Z(coefficients_for_Z),
            coefficients_for_H(coefficients_for_H) {}

    template<typename RingT>
    qrp_witness<RingT>::qrp_witness(const size_t num_variables,
                                    const size_t degree,
                                    const size_t num_inputs,
                                    const RingT &d1,
                                    const RingT &d2,
                                    const RingT &d3,
                                    const std::vector<RingT> &coefficients_for_ABCs,
                                    const std::vector<RingT> &coefficients_for_A_io,
                                    const std::vector<RingT> &coefficients_for_B_io,
                                    const std::vector<RingT> &coefficients_for_C_io,
                                    const std::vector<RingT> &coefficients_for_A_mid,
                                    const std::vector<RingT> &coefficients_for_B_mid,
                                    const std::vector<RingT> &coefficients_for_C_mid,
                                    std::vector<RingT> &&coefficients_for_H) :
            num_variables_(num_variables),
            degree_(degree),
            num_inputs_(num_inputs),
            d1(d1),
            d2(d2),
            d3(d3),
            coefficients_for_ABCs(coefficients_for_ABCs),
            coefficients_for_A_io(coefficients_for_A_io),
            coefficients_for_B_io(coefficients_for_B_io),
            coefficients_for_C_io(coefficients_for_C_io),
            coefficients_for_A_mid(coefficients_for_A_mid),
            coefficients_for_B_mid(coefficients_for_B_mid),
            coefficients_for_C_mid(coefficients_for_C_mid),
            coefficients_for_H(std::move(coefficients_for_H)) {}


    template<typename RingT>
    size_t qrp_witness<RingT>::num_variables() const {
        return num_variables_;
    }

    template<typename RingT>
    size_t qrp_witness<RingT>::degree() const {
        return degree_;
    }

    template<typename RingT>
    size_t qrp_witness<RingT>::num_inputs() const {
        return num_inputs_;
    }


} // ringsnark

#endif // QRP_TCC_
