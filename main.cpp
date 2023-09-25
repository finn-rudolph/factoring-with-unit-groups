#include <cassert>
#include <chrono>
#include <cstdint>
#include <future>
#include <iostream>
#include <variant>
#include <vector>

#include <gmpxx.h>

struct QuadFieldUnit
{
    mpz_class x, y;

    QuadFieldUnit() { x = 1, y = 0; }

    QuadFieldUnit(mpz_class &&x, mpz_class &&y) { this->x = x, this->y = y; }

    mpz_class norm(const mpz_class &d, const mpz_class &n)
    {
        return (x * x + d * y * y) % n;
    }

    static QuadFieldUnit square(const QuadFieldUnit &a, const mpz_class &d, const mpz_class &n)
    {
        return (QuadFieldUnit){((a.x * a.x) % n - (d * (a.y * a.y) % n) % n + n) % n,
                               ((a.x << 1) * a.y) % n};
    }

    static QuadFieldUnit mul(const QuadFieldUnit &a, const QuadFieldUnit &b, const mpz_class &d, const mpz_class &n)
    {
        mpz_class q = ((a.x + a.y) * (b.x + b.y)) % n,
                  r = (a.x * b.x) % n,
                  s = (a.y * b.y) % n;
        return QuadFieldUnit((r - (d * s) % n + n) % n, (q - r - s + (n << 1)) % n);
    }

    static QuadFieldUnit pow(QuadFieldUnit a, uint64_t exponent, const mpz_class &d, const mpz_class &n)
    {
        QuadFieldUnit result;
        while (exponent)
        {
            if (exponent & 1)
                result = mul(result, a, d, n);
            a = square(a, d, n);
            exponent >>= 1;
        }
        return result;
    }

    static std::variant<std::pair<QuadFieldUnit, mpz_class>, mpz_class> initial_point(
        const mpz_class &n, gmp_randclass &rng)
    {
        mpz_class x = rng.get_z_range(n), y = rng.get_z_range(n);

        mpz_class g, s, t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), y.get_mpz_t(), n.get_mpz_t());
        if (g > 1)
            return g;

        mpz_class d = ((((1 - x * x) % n) + n) * ((s * s) % n)) % n;
        return std::make_pair(QuadFieldUnit(std::move(x), std::move(y)), std::move(d));
    }

    void print()
    {
        std::cout << "Point2d { x = " << x << ", y = " << y << " }\n";
    }
};

// A unit modulo p in the cubic field Q(sqrt(d))
struct CubicFieldUnit
{
    mpz_class x, y, z;

    CubicFieldUnit() { x = 1, y = 0, z = 0; }

    CubicFieldUnit(mpz_class &&x, mpz_class &&y, mpz_class &&z)
    {
        this->x = x, this->y = y, this->z = z;
    }

    mpz_class norm(const mpz_class &d, const mpz_class &n)
    {
        return (((x * x) % n) * x +
                ((((y * y) % n) * y) % n + (((z * z) % n) * ((z * d) % n)) % n) * d -
                (3 * ((x * y) % n) * ((z * d) % n)) % n + n) %
               n;
    }

    static CubicFieldUnit mul(
        const CubicFieldUnit &a, const CubicFieldUnit &b, const mpz_class &d, const mpz_class &n)
    {
        const mpz_class s1 = (a.x * b.x) % n, s2 = (a.y * b.y) % n, s3 = (a.z * b.z) % n,
                        r1 = ((a.x + a.y) * (b.x + b.y)) % n, r2 = ((a.y + a.z) * (b.y + b.z)) % n,
                        r3 = ((a.z + a.x) * (b.z + b.x)) % n;
        return CubicFieldUnit((d * (r2 - s2 - s3 + 2 * n) + s1) % n,
                              (r1 - s1 - s2 + d * s3 + 2 * n) % n,
                              (r3 - s3 - s1 + s2 + 2 * n) % n);
    }

    static CubicFieldUnit square(const CubicFieldUnit &a, const mpz_class &d, const mpz_class &n)
    {
        return CubicFieldUnit((a.x * a.x + ((2 * d * a.y) % n) * a.z) % n,
                              (2 * a.x * a.y + d * ((a.z * a.z) % n)) % n,
                              (2 * a.z * a.x + a.y * a.y) % n);
    }

    static CubicFieldUnit pow(
        CubicFieldUnit a, uint64_t exponent, const mpz_class &d, const mpz_class &n)
    {
        CubicFieldUnit result;
        while (exponent)
        {
            if (exponent & 1)
                result = mul(result, a, d, n);
            a = square(a, d, n);
            exponent >>= 1;
        }
        return result;
    }

    static std::variant<std::pair<CubicFieldUnit, mpz_class>, mpz_class> initial_point(
        const mpz_class &n, gmp_randclass &rng)
    {
        mpz_class x = rng.get_z_range(n), y = rng.get_z_range(n);

        mpz_class g, s, t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), y.get_mpz_t(), n.get_mpz_t());
        if (g > 1)
            return g;

        const mpz_class d = ((1 - (((x * x) % n) * x) % n + n) * ((((s * s) % n) * (s % n + n)) % n)) % n;
        return std::make_pair(CubicFieldUnit(std::move(x), std::move(y), 0), std::move(d));
    }

    void print()
    {
        std::cout << "CubicFieldUnit { x = " << x << ", y = " << y << ", z = " << z << " }\n";
    }
};

struct Quaternion // <- misleading name, but it's close
{
    mpz_class x, y, z, w;

    Quaternion() { x = 1, y = 0, z = 0, w = 0; }

    Quaternion(mpz_class &&x, mpz_class &&y, mpz_class &&z, mpz_class &&w)
    {
        this->x = x, this->y = y, this->z = z, this->w = w;
    }

    mpz_class norm(const mpz_class &d, const mpz_class &n)
    {
        return (x * x + d * (y * y + z * z + w * w)) % n;
    }

    static Quaternion mul(const Quaternion &a, const Quaternion &b, const mpz_class &d, const mpz_class &n)
    {
        return Quaternion((a.x * b.x - ((d * (a.y * b.y + a.z * b.z + a.w * b.w)) % n) + n) % n,
                          (a.x * b.y + a.y * b.x + a.z * b.w - (a.w * b.z) % n + n) % n,
                          (a.x * b.z + a.z * b.x + a.w * b.y - (a.y * b.w) % n + n) % n,
                          (a.x * b.w + a.w * b.x + a.y * b.z - (a.z * b.y) % n + n) % n);
    }

    static Quaternion square(const Quaternion &a, const mpz_class &d, const mpz_class &n)
    {
        return Quaternion((a.x * a.x - ((d * (a.y * a.y + a.z * a.z + a.w * a.w)) % n) + n) % n,
                          ((a.x * a.y) << 1) % n, ((a.x * a.z) << 1) % n, ((a.x * a.w) << 1) % n);
    }

    static Quaternion pow(Quaternion a, uint64_t exponent, const mpz_class &d, const mpz_class &n)
    {
        Quaternion result;
        while (exponent)
        {
            if (exponent & 1)
                result = mul(result, a, d, n);
            a = square(a, d, n);
            exponent >>= 1;
        }
        return result;
    }

    static std::variant<std::pair<Quaternion, mpz_class>, mpz_class> initial_point(
        const mpz_class &n, gmp_randclass &rng)
    {
        mpz_class x = rng.get_z_range(n), y = rng.get_z_range(n), z = rng.get_z_range(n),
                  w = rng.get_z_range(n);

        mpz_class denom = (y * y + z * z + w * w) % n, g, s, t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), denom.get_mpz_t(), n.get_mpz_t());
        if (g > 1)
            return g;

        const mpz_class d = ((((1 - x * x) % n) + n) * (s % n + n)) % n;
        return std::make_pair(Quaternion(std::move(x), std::move(y), std::move(z), std::move(w)),
                              std::move(d));
    }

    void print()
    {
        std::cout << "Point4d { x = " << x << ", y = " << y
                  << ", z = " << z << ", w = " << w << " }\n";
    }
};

constexpr size_t SMOOTHNESS_BOUND_LIMIT = 1 << 24;
constexpr size_t NUM_TRIALS = 4;

template <typename T>
mpz_class factorize_with_algebraic_group(mpz_class const &n)
{
    std::vector<uint64_t> primes;

    {
        std::vector<bool> is_composite(SMOOTHNESS_BOUND_LIMIT + 1);
        for (size_t i = 2; i <= SMOOTHNESS_BOUND_LIMIT; ++i)
            if (!is_composite[i])
            {
                primes.push_back(i);
                for (size_t j = i * i; j <= SMOOTHNESS_BOUND_LIMIT; j += i)
                    is_composite[j] = 1;
            }
    }

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(std::chrono::duration_cast<std::chrono::nanoseconds>(
                 std::chrono::system_clock::now().time_since_epoch())
                 .count());

    for (size_t smoothness_bound = 2; smoothness_bound <= SMOOTHNESS_BOUND_LIMIT; smoothness_bound <<= 1)
    {
        for (size_t t = 0; t < NUM_TRIALS; ++t)
        {
            auto initial = T::initial_point(n, rng);
            if (std::holds_alternative<mpz_class>(initial))
                return std::get<mpz_class>(initial);

            auto [a, d] = std::get<std::pair<T, mpz_class>>(initial);
            assert(a.norm(d, n) == 1);

            for (uint64_t p : primes)
            {
                uint64_t q = p;
                while (q <= smoothness_bound)
                {
                    a = T::pow(a, p, d, n);
                    assert(a.norm(d, n) == 1);

                    mpz_class g;
                    mpz_gcd(g.get_mpz_t(), a.y.get_mpz_t(), n.get_mpz_t());
                    if (g > 1)
                        return g;
                    q *= p;
                }
            }
        }
    }

    return -1;
}

template <typename R>
bool is_ready(std::future<R> const &f)
{
    return f.wait_for(std::chrono::seconds(1)) == std::future_status::ready;
}

int main()
{
    std::cout << "n : ";
    mpz_class n;
    std::cin >> n;

    std::cout << "Do you want to work\n"
                 " - on the group of units in a quadratic field [p]\n"
                 " - on the group of units in a cubic field [c]\n"
                 " - or some quaternion-like thing [q]\n"
                 "=> [p / c / q] ";

    char method;
    std::cin >> method;

    std::vector<std::future<mpz_class>> fut;

    for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i)
    {
        if (method == 'p')
            fut.emplace_back(std::async(factorize_with_algebraic_group<QuadFieldUnit>, n));
        else if (method == 'c')
            fut.emplace_back(std::async(factorize_with_algebraic_group<CubicFieldUnit>, n));
        else
            fut.emplace_back(std::async(factorize_with_algebraic_group<Quaternion>, n));
    }

    while (!fut.empty())
    {
        for (auto f = fut.begin(); f != fut.end(); ++f)
            if (is_ready(*f))
            {
                mpz_class factor = f->get();
                if (factor != -1)
                {
                    std::cout << "factor found:\n"
                              << factor << '\n';
                    return 0;
                }
                else
                {
                    fut.erase(f);
                    break;
                }
            }
    }

    std::cout << "no factor found\n";
}