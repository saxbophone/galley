#ifndef COM_SAXBOPHONE_GALLEY_GF_HPP
#define COM_SAXBOPHONE_GALLEY_GF_HPP

#include <array>

#include <cstddef> // size_t
#include <cstdlib> // abort

namespace {
    // constexpr-safe partial reimplementation of std::optional<std::size_t>
    struct Nullable {
        constexpr Nullable() {}
        constexpr Nullable(std::size_t value) : _value(value), _present(true) {}
        constexpr Nullable& operator=(std::size_t value) {
            _value = value;
            _present = true;
            return *this;
        }
        constexpr std::size_t value() const {
            // lazy error-detection: abort() on request for missing value
            if (not _present) {
                abort();
            }
            return _value;
        }
        constexpr operator std::size_t() const {
            return value();
        }
        constexpr bool has_value() const {
            return _present;
        }
        constexpr operator bool() const {
            return has_value();
        }
        constexpr void reset() {
            _value = 0;
            _present = false;
        }
    private:
        std::size_t _value = 0;
        bool _present = false;
    };

    // NOTE: N must be prime!
    // currently not enforced, but should be enforced later, ideally at template-level
    template <std::size_t N>
    class GaloisFieldLookupTable {
    public:
        // TODO: generate compile-time galois field lookup table here!
        constexpr GaloisFieldLookupTable() {
            for (std::size_t i = 1; i < N; i++) {
                // additive inverse
                _additive_inverse[i] = N - i;
                // square
                _square[i] = (i * i) % N;
                // multiplicative inverse --use the "find the missing term" trick
                // x = 1/i mod N
                // x = (y*N + 1) / i
                for (std::size_t y = 0; y < N; y++) {
                    std::size_t candidate = (y * N + 1) / i;
                    // test candidate before accepting
                    if ((i * candidate) % N == 1) {
                        _multiplicative_inverse[i] = candidate;
                        break;
                    }
                }
            }
            for (std::size_t i = 1; i < N; i++) {
                // square root can be populated by the reverse of square. just do linear search for it
                for (std::size_t j = 1; j < N; j++) {
                    if (_square[j].value() == i) {
                        _square_root[i] = j;
                        break;
                    }
                }
            }
            // additive inverse zero is a special case, we also do square in this way
            _additive_inverse[0] = 0;
            _square[0] = 0;
        }
        constexpr Nullable additive_inverse(std::size_t x) const {
            return _additive_inverse[x];
        }
        constexpr Nullable multiplicative_inverse(std::size_t x) const {
            return _multiplicative_inverse[x];
        }
        constexpr Nullable square(std::size_t x) const {
            return _square[x];
        }
        constexpr Nullable square_root(std::size_t x) const {
            return _square_root[x];
        }
    private:
        std::array<Nullable, N> _additive_inverse;
        std::array<Nullable, N> _multiplicative_inverse;
        std::array<Nullable, N> _square;
        std::array<Nullable, N> _square_root;
    };
}

namespace com::saxbophone::galley {
    // NOTE: N must be prime!
    // currently not enforced, but should be enforced later, ideally at template-level
    template <std::size_t N>
    class GF {
    public:
        constexpr GF() {}
        constexpr GF(std::size_t x) : _value(x % N) {}
        constexpr std::size_t get_value() const {
            return _value;
        }
        constexpr operator std::size_t() const {
            return get_value();
        }
        constexpr void set_value(std::size_t x) {
            // NOTE: this should always be reduced mod-N before being stored!
            // XXX: out-of-bounds values are wrapped silently --do we want that?
            _value = x % N;
        }
        constexpr GF& operator=(std::size_t x) {
            set_value(x);
            return *this;
        }
        constexpr GF& operator+=(GF rhs) {
            set_value(_value + rhs._value);
            return *this;
        }
        friend constexpr GF operator+(GF lhs, const GF& rhs) {
            lhs += rhs; // reuse compound assignment
            return lhs; // return the result by value (uses move constructor)
        }
        // unary minus operator, aka finite field additive inverse
        constexpr GF operator-() const {
            return _LOOKUP_TABLE.additive_inverse(_value).value();
        }
        constexpr GF& operator-=(GF rhs) {
            // here's where we apply some finite field arithmetic magic!
            set_value(_value + (-rhs)._value);
            return *this;
        }
        friend constexpr GF operator-(GF lhs, const GF& rhs) {
            lhs -= rhs; // reuse compound assignment
            return lhs; // return the result by value (uses move constructor)
        }
        constexpr GF& operator*=(GF rhs) {
            set_value(_value * rhs._value);
            return *this;
        }
        friend constexpr GF operator*(GF lhs, const GF& rhs) {
            lhs *= rhs; // reuse compound assignment
            return lhs; // return the result by value (uses move constructor)
        }
        constexpr GF& operator/=(GF rhs) {
            // here's where we apply some finite field arithmetic magic!
            set_value(*this * inv(rhs));
            return *this;
        }
        friend constexpr GF operator/(GF lhs, const GF& rhs) {
            lhs /= rhs; // reuse compound assignment
            return lhs; // return the result by value (uses move constructor)
        }
        static constexpr GF squared(const GF& x) {
            return _LOOKUP_TABLE.square(x._value).value();
        }
        static constexpr GF sqrt(const GF& x) {
            return _LOOKUP_TABLE.square_root(x._value).value();
        }
        // NOTE: just exposes multiplicative inverse (x⁻¹) in case it's wanted directly
        static constexpr GF inv(const GF& x) {
            return _LOOKUP_TABLE.multiplicative_inverse(x._value).value();
        }
    private:
        // "= {}" in the definition should be redundant, but GCC has a bug where it doesn't abide by P0386, requiring this
        static constexpr GaloisFieldLookupTable<N> _LOOKUP_TABLE = {};
        std::size_t _value = 0;
    };
} // namespace com::saxbophone::galley
#endif // include guard
