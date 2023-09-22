#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <variant>
#include <vector>

#include <gmpxx.h>

struct Point2d
{
    mpz_class x, y;

    Point2d() { x = 1, y = 0; }

    Point2d(mpz_class &&x, mpz_class &&y) { this->x = x, this->y = y; }

    mpz_class norm(const mpz_class &d, const mpz_class &n)
    {
        return (x * x + d * y * y) % n;
    }

    static Point2d square(const Point2d &a, const mpz_class &d, const mpz_class &n)
    {
        return (Point2d){((a.x * a.x) % n - (d * (a.y * a.y) % n) % n + n) % n,
                         ((a.x << 1) * a.y) % n};
    }

    static Point2d mul(const Point2d &a, const Point2d &b, const mpz_class &d, const mpz_class &n)
    {
        mpz_class q = ((a.x + a.y) * (b.x + b.y)) % n,
                  r = (a.x * b.x) % n,
                  s = (a.y * b.y) % n;
        return Point2d((r - (d * s) % n + n) % n, (q - r - s + (n << 1)) % n);
    }

    static Point2d pow(Point2d a, uint64_t exponent, const mpz_class &d, const mpz_class &n)
    {
        Point2d result;
        while (exponent)
        {
            if (exponent & 1)
                result = mul(result, a, d, n);
            a = square(a, d, n);
            exponent >>= 1;
        }
        return result;
    }

    static std::variant<std::pair<Point2d, mpz_class>, mpz_class> initial_point(
        const mpz_class &n, gmp_randclass &rng)
    {
        mpz_class x = rng.get_z_range(n), y = rng.get_z_range(n);

        mpz_class g, s, t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), y.get_mpz_t(), n.get_mpz_t());
        if (g > 1)
            return g;

        mpz_class d = ((((1 - x * x) % n) + n) * ((s * s) % n)) % n;
        return std::make_pair(Point2d(std::move(x), std::move(y)), std::move(d));
    }

    void print()
    {
        std::cout << "Point2d { x = " << x << ", y = " << y << " }\n";
    }
};

struct Point4d
{
    mpz_class x, y, z, w;

    Point4d() { x = 1, y = 0, z = 0, w = 0; }

    Point4d(mpz_class &&x, mpz_class &&y, mpz_class &&z, mpz_class &&w)
    {
        this->x = x, this->y = y, this->z = z, this->w = w;
    }

    mpz_class norm(const mpz_class &d, const mpz_class &n)
    {
        return (x * x + d * (y * y + z * z + w * w)) % n;
    }

    static Point4d mul(const Point4d &a, const Point4d &b, const mpz_class &d, const mpz_class &n)
    {
        return Point4d((a.x * b.x - ((d * (a.y * b.y + a.z * b.z + a.w * b.w)) % n) + n) % n,
                       (a.x * b.y + a.y * b.x + a.z * b.w - (a.w * b.z) % n + n) % n,
                       (a.x * b.z + a.z * b.x + a.w * b.y - (a.y * b.w) % n + n) % n,
                       (a.x * b.w + a.w * b.x + a.y * b.z - (a.z * b.y) % n + n) % n);
    }

    static Point4d square(const Point4d &a, const mpz_class &d, const mpz_class &n)
    {
        return Point4d((a.x * a.x - ((d * (a.y * a.y + a.z * a.z + a.w * a.w)) % n) + n) % n,
                       ((a.x * a.y) << 1) % n, ((a.x * a.z) << 1) % n, ((a.x * a.w) << 1) % n);
    }

    static Point4d pow(Point4d a, uint64_t exponent, const mpz_class &d, const mpz_class &n)
    {
        Point4d result;
        while (exponent)
        {
            if (exponent & 1)
                result = mul(result, a, d, n);
            a = square(a, d, n);
            exponent >>= 1;
        }
        return result;
    }

    static std::variant<std::pair<Point4d, mpz_class>, mpz_class> initial_point(
        const mpz_class &n, gmp_randclass &rng)
    {
        mpz_class x = rng.get_z_range(n), y = rng.get_z_range(n), z = rng.get_z_range(n),
                  w = rng.get_z_range(n);

        mpz_class denom = (y * y + z * z + w * w) % n, g, s, t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), denom.get_mpz_t(), n.get_mpz_t());
        if (g > 1)
            return g;

        const mpz_class d = ((((1 - x * x) % n) + n) * (s % n + n)) % n;
        return std::make_pair(Point4d(std::move(x), std::move(y), std::move(z), std::move(w)),
                              std::move(d));
    }

    void print()
    {
        std::cout << "Point4d { x = " << x << ", y = " << y
                  << ", z = " << z << ", w = " << w << " }\n";
    }
};

constexpr size_t SMOOTHNESS_BOUND_LIMIT = 1 << 20;
constexpr size_t NUM_TRIALS = 32;

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
    rng.seed(std::chrono::duration_cast<std::chrono::microseconds>(
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

int main()
{
    std::cout << "n : ";
    mpz_class n;
    std::cin >> n;

    std::cout << "Do you want to work on the group of units in a quadratic field [p]"
                 "or some quaternion-like thing [q]? [p / q] ";

    char method;
    std::cin >> method;

    mpz_class factor;
    if (method == 'p')
        factor = factorize_with_algebraic_group<Point2d>(n);
    else
        factor = factorize_with_algebraic_group<Point4d>(n);

    if (factor == -1)
        std::cout << "no factor found\n";
    else
        std::cout << "factor found:\n"
                  << factor << '\n';
}