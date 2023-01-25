#ifndef RINGSNARK_TEST_UTILS_HPP
#define RINGSNARK_TEST_UTILS_HPP

#include <random>

namespace {
    template<typename T>
    class PrimitiveWrapper {
    public:
        T val;

        static PrimitiveWrapper zero() { return PrimitiveWrapper(0); }

        static PrimitiveWrapper one() { return PrimitiveWrapper(1); }

        explicit PrimitiveWrapper(T val = 0) : val(val) {}

        bool is_zero() const { return val == 0; }

        inline static std::default_random_engine e = std::default_random_engine(42);

        static PrimitiveWrapper<T> random_element() {
            using namespace std;
            uniform_int_distribution<size_t> u(0, 1000);
            size_t rand = u(e);
            return PrimitiveWrapper<T>((double) rand);
        }

        static PrimitiveWrapper<T> random_exceptional_element() {
            using namespace std;
            uniform_int_distribution<size_t> u(0, 10);
            size_t rand = u(e);
            return PrimitiveWrapper<T>((double) rand);
        }

        PrimitiveWrapper &operator-() const {
            return *(new PrimitiveWrapper(-this->val));
        }

        PrimitiveWrapper &operator+=(const PrimitiveWrapper &rhs) {
            val += rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator-=(const PrimitiveWrapper &rhs) {
            val -= rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator*=(const PrimitiveWrapper &rhs) {
            val *= rhs.val;
            return *this;
        }

        PrimitiveWrapper &operator/=(const PrimitiveWrapper &rhs) {
            val /= rhs.val;
            return *this;
        }

        bool operator==(const PrimitiveWrapper<T> &rhs) const {
            return abs(val - rhs.val) <= 1e-3 * abs(val); // Approximate equality for double
        }

        bool operator!=(const PrimitiveWrapper<T> &rhs) const {
            return !this->operator==(rhs);
        }
    };

    template<typename T>
    inline PrimitiveWrapper<T> operator+(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator-(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator*(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    inline PrimitiveWrapper<T> operator/(PrimitiveWrapper<T> lhs, const PrimitiveWrapper<T> &rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T>
    inline std::ostream &operator<<(std::ostream &out, PrimitiveWrapper<T> const &w) {
        return out << w.val;
    }
}

#endif //RINGSNARK_TEST_UTILS_HPP
