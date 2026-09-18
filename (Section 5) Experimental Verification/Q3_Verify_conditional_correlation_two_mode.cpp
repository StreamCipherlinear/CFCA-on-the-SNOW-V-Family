/**
 * verify_correlation.cpp
 *
 * Compares theoretical conditional correlation with empirical measurement
 * using theoretical standard deviation from CLT.
 *
 * Usage: ./verify_correlation [--seed N] [--symbols N] [--pair-symbols N]
 *                             [--repetitions N] [--B N] [--L N] [--K N] [--ps P]
 *   --seed N         : fixed seed (random if absent)
 *   --symbols N      : symbols in direct mode       (default: theoretical M)
 *   --pair-symbols N : records in collision mode    (default: sqrt(M_pair*2^(K+1)))
 *   --repetitions N  : number of independent runs   (default 100)
 *   --B N            : unknown bits, direct mode    (default 16)
 *   --L N            : unknown bits, collision mode (default 24)
 *   --K N            : bits cancelled by collisions (default 16)
 *   --ps P           : target success probability   (default 0.999992)
 *
 * M follows Q4: M = (Phi^{-1}(ps) + sqrt(2 ln2 B_eff))^2 / rho^2, with
 * rho^2 = aver_rho_sq in direct mode and aver_rho_sq^2 in collision mode.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <tuple>
#include <vector>
#include <random>
#include <unordered_map>

// ----------------------------------------------------------------------
// Bit manipulation helpers
// ----------------------------------------------------------------------

static inline unsigned parity64(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_parityll(value));
#else
    value ^= value >> 32; value ^= value >> 16;
    value ^= value >> 8;  value ^= value >> 4;
    value &= 0x0FU;
    return static_cast<unsigned>((0x6996U >> value) & 1U);
#endif
}

static inline unsigned dot16(uint16_t mask, uint16_t value) {
    return parity64(static_cast<uint64_t>(mask & value));
}

static inline uint8_t gf16_mul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (unsigned i = 0; i < 4; ++i) {
        if (b & 1U) p ^= a;
        a <<= 1;
        if (a & 0x10U) a ^= 0x13U;
        b >>= 1;
    }
    return p & 0x0FU;
}

static const uint8_t SBOX[16] = {
        0x9,0x4,0xA,0xB, 0xD,0x1,0x8,0x5,
        0x6,0x2,0x0,0x3, 0xC,0xE,0xF,0x7
};

static inline uint8_t sr(uint8_t v) { return SBOX[v & 0x0FU]; }

static uint8_t S8_transform(uint8_t value) {
    uint8_t x0 = value & 0x0FU, x1 = (value >> 4) & 0x0FU;
    uint8_t y0 = sr(x0), y1 = sr(x1);
    uint8_t z0 = y0 ^ gf16_mul(4, y1);
    uint8_t z1 = gf16_mul(4, y0) ^ y1;
    return (z0 & 0x0FU) | ((z1 & 0x0FU) << 4);
}

static inline uint16_t S16_transform(uint16_t value) {
    uint8_t hi = S8_transform(static_cast<uint8_t>(value >> 8));
    uint8_t lo = S8_transform(static_cast<uint8_t>(value));
    return (static_cast<uint16_t>(hi) << 8) | lo;
}

// ----------------------------------------------------------------------
// LFSR and symbolic representation
// ----------------------------------------------------------------------

static const unsigned T1_SHIFT = 32;
static const unsigned T2_SHIFT = 0;
static const uint64_t LFSR_REDUCTION = UINT64_C(0x1B);

struct Lfsr { uint32_t high; uint32_t low; };

static inline uint64_t get_state(const Lfsr &l) {
    return (static_cast<uint64_t>(l.high) << 32) | l.low;
}

static inline void set_state(Lfsr &l, uint64_t s) {
    if (s == 0) s = 1;
    l.high = static_cast<uint32_t>(s >> 32);
    l.low  = static_cast<uint32_t>(s);
}

static inline uint16_t get_tap(const Lfsr &l, unsigned sh) {
    return static_cast<uint16_t>((get_state(l) >> sh) & 0xFFFFU);
}

static inline uint16_t get_T1(const Lfsr &l) { return get_tap(l, T1_SHIFT); }
static inline uint16_t get_T2(const Lfsr &l) { return get_tap(l, T2_SHIFT); }

static inline uint64_t lfsr_clock_bit_value(uint64_t s) {
    uint64_t out = s >> 63;
    s <<= 1;
    if (out) s ^= LFSR_REDUCTION;
    return s;
}

static inline void lfsr_clock_bit(Lfsr &l) {
    set_state(l, lfsr_clock_bit_value(get_state(l)));
}

static inline void lfsr_clock_word(Lfsr &l) {
    for (unsigned i = 0; i < 16; ++i) lfsr_clock_bit(l);
}

struct SymbolicLfsr { std::array<uint64_t, 64> row; };

static SymbolicLfsr make_symbolic_lfsr() {
    SymbolicLfsr s;
    for (unsigned i = 0; i < 64; ++i) s.row[i] = UINT64_C(1) << i;
    return s;
}

static inline void symbolic_clock_bit(SymbolicLfsr &s) {
    uint64_t out = s.row[63];
    for (int i = 63; i >= 1; --i) s.row[i] = s.row[i-1];
    s.row[0] = out;
    s.row[1] ^= out; s.row[3] ^= out; s.row[4] ^= out;
}

static inline void symbolic_clock_word(SymbolicLfsr &s) {
    for (unsigned i = 0; i < 16; ++i) symbolic_clock_bit(s);
}

static uint64_t apply_symbolic_mask(const std::array<uint64_t, 16> &masks, uint16_t mask) {
    uint64_t r = 0;
    for (unsigned i = 0; i < 16; ++i) if ((mask >> i) & 1U) r ^= masks[i];
    return r;
}

// ----------------------------------------------------------------------
// Cipher state and sample generation
// ----------------------------------------------------------------------

struct Cipher {
    Lfsr lfsr;
    SymbolicLfsr symbolic;
    uint16_t R1, R2, R3;
};

struct Sample {
    uint16_t T1, T2, z;
    std::array<uint64_t, 16> T1_bit_masks, T2_bit_masks;
};

static inline void fsm_update(Cipher &c, uint16_t T2) {
    uint16_t o1 = c.R1, o2 = c.R2, o3 = c.R3;
    c.R1 = static_cast<uint16_t>(o2 + static_cast<uint16_t>(o3 ^ T2));
    c.R2 = S16_transform(o1);
    c.R3 = S16_transform(o2);
}

static Sample generate_sample(Cipher &c) {
    Sample s;
    s.T1 = get_T1(c.lfsr);
    s.T2 = get_T2(c.lfsr);
    s.z  = static_cast<uint16_t>(static_cast<uint16_t>(c.R1 + s.T1) ^ c.R2);
    for (unsigned i = 0; i < 16; ++i) {
        s.T1_bit_masks[i] = c.symbolic.row[T1_SHIFT + i];
        s.T2_bit_masks[i] = c.symbolic.row[T2_SHIFT + i];
    }
    fsm_update(c, s.T2);
    lfsr_clock_word(c.lfsr);
    symbolic_clock_word(c.symbolic);
    return s;
}

static uint64_t splitmix64(uint64_t &state) {
    uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

static Cipher initialize_cipher(uint64_t seed) {
    Cipher c;
    uint64_t st = splitmix64(seed);
    if (st == 0) st = 1;
    set_state(c.lfsr, st);
    c.R1 = static_cast<uint16_t>(splitmix64(seed));
    c.R2 = static_cast<uint16_t>(splitmix64(seed));
    c.R3 = static_cast<uint16_t>(splitmix64(seed));
    c.symbolic = make_symbolic_lfsr();
    return c;
}

// ----------------------------------------------------------------------
// Automaton8 correlation function
// ----------------------------------------------------------------------
static double automaton8(int input1, int input2, int output) {
    double k = 0.0;
    int e = 0;
    double sign = 1.0;
    for (int t = 0; t < 8; ++t) {
        int flag = ((output >> (7 - t)) & 1) * 4 +
                   ((input1 >> (7 - t)) & 1) * 2 +
                   ((input2 >> (7 - t)) & 1);
        if (e == 0) {
            if (flag == 7)      e = 1;
            else if (flag == 0) e = 0;
            else                return -100.0;
        } else {
            if (flag == 0 || flag == 3 || flag == 5 || flag == 6) {
                e = 1; k += 1.0;
            } else {
                e = 0; k += 1.0;
            }
            if (flag == 3 || flag == 4) sign = -sign;
        }
    }
    return sign * std::pow(2.0, -k);
}

// ----------------------------------------------------------------------
// Adaptive approximation table
// ----------------------------------------------------------------------

struct AdaptiveEntry { uint16_t l; unsigned tau; double absolute_correlation; };

static void fwht_small(int values[256]) {
    for (unsigned len = 1; len < 256; len <<= 1)
        for (unsigned st = 0; st < 256; st += len << 1)
            for (unsigned o = 0; o < len; ++o) {
                int a = values[st + o], b = values[st + o + len];
                values[st + o] = a + b;
                values[st + o + len] = a - b;
            }
}

static std::array<AdaptiveEntry, 16> compute_adaptive_table(
        const std::vector<int>& xi_list,
        uint16_t beta,
        uint16_t gamma,
        uint16_t h,
        uint16_t n)
{
    std::array<AdaptiveEntry, 16> table;
    const double corr_e1_abs = 0.25;

    std::vector<double> C(xi_list.size());
    for (size_t idx = 0; idx < xi_list.size(); ++idx) {
        int xi = xi_list[idx];
        double sum = 0.0;
        for (int theta = 0; theta < 256; ++theta) {
            double v1 = automaton8(theta, static_cast<int>(n & 0xFF), static_cast<int>(gamma & 0xFF));
            double v2 = automaton8((beta & 0xFF) ^ xi, static_cast<int>(h & 0xFF), theta);
            if (v1 > -100.0 && v2 > -100.0)
                sum += v1 * v2;
        }
        C[idx] = sum;
    }

    for (unsigned u = 0; u < 16; ++u) {
        double combined[256] = {0.0};
        for (size_t idx = 0; idx < xi_list.size(); ++idx) {
            int xi = xi_list[idx];
            int approx_x[256], approx_z[256];
            for (int v = 0; v < 256; ++v) {
                approx_x[v] = dot16(static_cast<uint16_t>(xi), S16_transform(static_cast<uint16_t>(v)));
                approx_z[v] = dot16(static_cast<uint16_t>(h & 0xFF), S16_transform(static_cast<uint16_t>(v)));
            }
            int spectrum[256];
            for (int y = 0; y < 256; ++y) {
                int score = 0;
                for (int x = 0; x < 256; ++x) {
                    int z = u ^ ((x + y) & 0xFF);
                    int bit = approx_x[x] ^ approx_z[z];
                    score += bit ? -1 : 1;
                }
                spectrum[y] = score;
            }
            fwht_small(spectrum);
            for (int l = 0; l < 256; ++l) {
                double rho = static_cast<double>(spectrum[l]) / 65536.0;
                combined[l] += rho * C[idx];
            }
        }
        int best_l = 0;
        double best_score = combined[0];
        for (int l = 1; l < 256; ++l) {
            if (std::fabs(combined[l]) > std::fabs(best_score)) {
                best_score = combined[l];
                best_l = l;
            }
        }
        AdaptiveEntry e;
        e.l = static_cast<uint16_t>(best_l);
        e.tau = (best_score < 0.0) ? 0U : 1U;
        e.absolute_correlation = corr_e1_abs * std::fabs(best_score);
        table[u] = e;
    }
    return table;
}

// ----------------------------------------------------------------------
// Helper: run one experiment
// ----------------------------------------------------------------------
static std::tuple<double, double, uint64_t> run_one_experiment(
        const std::array<AdaptiveEntry, 16>& table,
        uint64_t seed,
        uint64_t num_symbols,
        uint64_t num_symbols_pair,
        const uint16_t beta,
        const uint16_t gamma,
        const uint16_t m,
        const uint16_t n,
        const uint16_t h)
{
    // Direct mode
    Cipher cipher_dir = initialize_cipher(seed);
    uint64_t init_state = get_state(cipher_dir.lfsr);
    Sample prev, cur, next;
    prev = generate_sample(cipher_dir);
    cur  = generate_sample(cipher_dir);
    next = generate_sample(cipher_dir);

    uint64_t total_noise_dir = 0;
    uint64_t total_eq_dir = 0;

    for (uint64_t idx = 0; idx < num_symbols; ++idx) {
        unsigned u = prev.z & 0x0F;
        const AdaptiveEntry &e = table[u];

        unsigned observed = dot16(beta, cur.z) ^ dot16(gamma, next.z) ^ e.tau;
        uint64_t state_mask = apply_symbolic_mask(prev.T1_bit_masks, e.l) ^
                              apply_symbolic_mask(cur.T1_bit_masks, m) ^
                              apply_symbolic_mask(next.T1_bit_masks, n) ^
                              apply_symbolic_mask(cur.T2_bit_masks, h);

        unsigned noise = observed ^ parity64(state_mask & init_state);
        total_noise_dir += noise;
        total_eq_dir += 1;

        prev = cur;
        cur = next;
        next = generate_sample(cipher_dir);
    }

    double p_dir = (double)total_noise_dir / total_eq_dir;
    double emp_dir_abs = std::fabs(1.0 - 2.0 * p_dir);

    // Pair mode
    uint64_t collision_mask = 0xFFFFULL << 48;
    std::unordered_map<uint64_t, std::vector<std::pair<uint64_t, unsigned>>> buckets;

    Cipher cipher_pair = initialize_cipher(seed + 1);
    uint64_t init_state_pair = get_state(cipher_pair.lfsr);

    Sample pv, cu, nx;
    pv = generate_sample(cipher_pair);
    cu = generate_sample(cipher_pair);
    nx = generate_sample(cipher_pair);

    uint64_t collected = 0;
    while (collected < num_symbols_pair) {
        unsigned u = pv.z & 0x0F;
        const AdaptiveEntry &e = table[u];
        unsigned observed = dot16(beta, cu.z) ^ dot16(gamma, nx.z) ^ e.tau;
        uint64_t state_mask = apply_symbolic_mask(pv.T1_bit_masks, e.l) ^
                              apply_symbolic_mask(cu.T1_bit_masks, m) ^
                              apply_symbolic_mask(nx.T1_bit_masks, n) ^
                              apply_symbolic_mask(cu.T2_bit_masks, h);

        uint64_t key = (state_mask & collision_mask) >> 48;
        buckets[key].push_back({state_mask, observed});
        collected++;

        pv = cu; cu = nx; nx = generate_sample(cipher_pair);
    }

    uint64_t total_noise_pair = 0;
    uint64_t total_pairs = 0;
    for (auto &kv : buckets) {
        auto &vec = kv.second;
        size_t sz = vec.size();
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = i + 1; j < sz; ++j) {
                const auto &a = vec[i];
                const auto &b = vec[j];
                uint64_t delta_mask = a.first ^ b.first;
                unsigned delta_obs = a.second ^ b.second;
                unsigned noise = delta_obs ^ parity64(delta_mask & init_state_pair);
                total_noise_pair += noise;
                total_pairs++;
            }
        }
    }

    double emp_pair_abs = 0.0;
    if (total_pairs > 0) {
        double p_pair = (double)total_noise_pair / total_pairs;
        emp_pair_abs = std::fabs(1.0 - 2.0 * p_pair);
    }

    return {emp_dir_abs, emp_pair_abs, total_pairs};
}

// ----------------------------------------------------------------------
// Normal CDF inverse (from Q4)
// ----------------------------------------------------------------------
static double inverse_normal_cdf(double p) {
    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return  std::numeric_limits<double>::infinity();
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                               -2.759285104469687e+02, 1.383577518672690e+02,
                               -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                               -1.556989798598866e+02, 6.680131188771972e+02,
                               -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03,-3.223964580411365e-01,
                               -2.400758277161838e+00,-2.549732539343734e+00,
                                4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = { 7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};
    double low = 0.02425, high = 1.0 - low;
    if (p < low) {
        double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
    if (p > high) {
        double q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0));
    }
    double q = p - 0.5, r = q * q;
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
}

// ----------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------
int main(int argc, char **argv) {
    uint64_t seed = 0;
    bool seed_provided = false;
    uint64_t symbols_direct = 0;       // 0 = auto
    uint64_t symbols_pair = 0;         // 0 = auto
    unsigned recover_bits = 16;        // direct mode: unknown bits
    unsigned recover_bits_pair = 24;   // collision mode: unknown bits
    unsigned collision_bits = 16;      // bits cancelled by collisions
    double success_probability = 0.999992;
    unsigned repetitions = 100;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--seed") == 0 && i+1 < argc) {
            seed = strtoull(argv[++i], NULL, 0);
            seed_provided = true;
        } else if (strcmp(argv[i], "--symbols") == 0 && i+1 < argc) {
            symbols_direct = strtoull(argv[++i], NULL, 0);
        } else if (strcmp(argv[i], "--pair-symbols") == 0 && i+1 < argc) {
            symbols_pair = strtoull(argv[++i], NULL, 0);
        } else if (strcmp(argv[i], "--B") == 0 && i+1 < argc) {
            recover_bits = static_cast<unsigned>(strtoul(argv[++i], NULL, 0));
        } else if (strcmp(argv[i], "--L") == 0 && i+1 < argc) {
            recover_bits_pair = static_cast<unsigned>(strtoul(argv[++i], NULL, 0));
        } else if (strcmp(argv[i], "--K") == 0 && i+1 < argc) {
            collision_bits = static_cast<unsigned>(strtoul(argv[++i], NULL, 0));
        } else if (strcmp(argv[i], "--ps") == 0 && i+1 < argc) {
            success_probability = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--repetitions") == 0 && i+1 < argc) {
            repetitions = static_cast<unsigned>(strtoul(argv[++i], NULL, 0));
            if (repetitions == 0) repetitions = 1;
        } else if (strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s [--seed N] [--symbols N] [--pair-symbols N] [--repetitions N]"
                   " [--B N] [--L N] [--K N] [--ps P]\n", argv[0]);
            printf("  --symbols N      : symbols in direct mode (default: theoretical M)\n");
            printf("  --pair-symbols N : records in collision mode (default: sqrt(M*2^(K+1)))\n");
            printf("  --B N / --L N / --K N / --ps P : parameters of the M expression\n");
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    // Fixed parameters
    const uint16_t beta_fixed = 0x0008;
    const uint16_t gamma_fixed = 0x0041;
    const uint16_t h_fixed = 0x0041;
    const uint16_t n_fixed = 0x0061;
    const uint16_t m_fixed = 0x0008;

    std::vector<int> xi_list = {
            0x17, 0x28, 0x3F, 0x41, 0x56, 0x69, 0x7E,
            0x83, 0x94, 0xAB, 0xBC, 0xC2, 0xD5, 0xEA, 0xFD
    };

    auto table = compute_adaptive_table(xi_list, beta_fixed, gamma_fixed, h_fixed, n_fixed);

    // Theoretical values
    double aver_rho_sq = 0.0;
    for (auto &e : table) aver_rho_sq += e.absolute_correlation * e.absolute_correlation;
    aver_rho_sq /= 16.0;
    double rms_direct = std::sqrt(aver_rho_sq);

    // Theoretical M, Q4 expression (B, not B+1)
    double z = inverse_normal_cdf(success_probability);
    double sqrt_term_direct = std::sqrt(2.0 * std::log(2.0) * recover_bits);
    long double M_direct = (1.0L / (long double)aver_rho_sq) *
                           std::pow((long double)(z + sqrt_term_direct), 2);

    unsigned eff_pair = (recover_bits_pair > collision_bits)
                        ? (recover_bits_pair - collision_bits) : 0u;
    double sqrt_term_pair = std::sqrt(2.0 * std::log(2.0) * eff_pair);
    long double M_pair = (1.0L / (long double)(aver_rho_sq * aver_rho_sq)) *
                         std::pow((long double)(z + sqrt_term_pair), 2);

    uint64_t num_symbols = (symbols_direct != 0)
                           ? symbols_direct
                           : static_cast<uint64_t>(std::ceil(M_direct));
    uint64_t num_symbols_pair = (symbols_pair != 0)
                                ? symbols_pair
                                : static_cast<uint64_t>(std::ceil(
                                      std::sqrt((long double)M_pair *
                                                std::ldexp(1.0L, (int)collision_bits + 1))));

    printf("===== Theoretical vs Empirical Correlation (with Theoretical Sigma) =====\n");
    printf("ps = %.6f, repetitions = %u\n", success_probability, repetitions);
    printf("aver_rho_sq = %.9f\n", aver_rho_sq);
    printf("Theoretical RMS |cor| (direct) = %.9f\n", rms_direct);
    printf("Direct    : B = %u, M = %.2f (2^%.4f)\n",
           recover_bits, (double)M_direct, (double)std::log2((double)M_direct));
    printf("Collision : L = %u, K = %u, M = %.2f (2^%.4f), records = %llu (2^%.4f)\n",
           recover_bits_pair, collision_bits, (double)M_pair,
           (double)std::log2((double)M_pair),
           (unsigned long long)num_symbols_pair,
           (double)std::log2((double)num_symbols_pair));
    if (symbols_direct != 0)
        printf("  (direct symbols overridden to %llu)\n", (unsigned long long)num_symbols);
    if (symbols_pair != 0)
        printf("  (collision records overridden to %llu)\n", (unsigned long long)num_symbols_pair);

    double sigma_dir = 1.0 / std::sqrt((double)num_symbols);
    double two_sigma_dir = 2.0 * sigma_dir;
    double three_sigma_dir = 3.0 * sigma_dir;

    int pass_2_dir = 0, pass_3_dir = 0;
    int pass_2_pair = 0, pass_3_pair = 0;

    std::vector<double> emp_dir_list, emp_pair_list;
    emp_dir_list.reserve(repetitions);
    emp_pair_list.reserve(repetitions);

    for (unsigned rep = 0; rep < repetitions; ++rep) {
        uint64_t current_seed;
        if (seed_provided) {
            current_seed = seed + rep;
        } else {
            std::random_device rd;
            current_seed = (static_cast<uint64_t>(rd()) << 32) | rd();
            if (current_seed == 0) current_seed = 1;
        }

        auto [emp_dir, emp_pair, pair_count] = run_one_experiment(
                table, current_seed, num_symbols, num_symbols_pair,
                beta_fixed, gamma_fixed, m_fixed, n_fixed, h_fixed);

        emp_dir_list.push_back(emp_dir);
        emp_pair_list.push_back(emp_pair);

        double diff_dir = emp_dir - rms_direct;
        if (std::fabs(diff_dir) <= 2.0 * sigma_dir) pass_2_dir++;
        if (std::fabs(diff_dir) <= 3.0 * sigma_dir) pass_3_dir++;

        double sigma_pair = 1.0 / std::sqrt((double)pair_count);
        double diff_pair = emp_pair - aver_rho_sq;
        if (std::fabs(diff_pair) <= 2.0 * sigma_pair) pass_2_pair++;
        if (std::fabs(diff_pair) <= 3.0 * sigma_pair) pass_3_pair++;
    }

    double min_dir = *std::min_element(emp_dir_list.begin(), emp_dir_list.end());
    double max_dir = *std::max_element(emp_dir_list.begin(), emp_dir_list.end());
    double min_pair = *std::min_element(emp_pair_list.begin(), emp_pair_list.end());
    double max_pair = *std::max_element(emp_pair_list.begin(), emp_pair_list.end());

    double mean_dir = 0.0, mean_pair = 0.0;
    for (size_t i = 0; i < emp_dir_list.size(); ++i) {
        mean_dir += emp_dir_list[i];
        mean_pair += emp_pair_list[i];
    }
    mean_dir /= repetitions;
    mean_pair /= repetitions;

    printf("\n===== Summary over %u repetitions =====\n", repetitions);
    printf("\nDirect mode:\n");
    printf("  Theoretical |cor|       : %.9f\n", rms_direct);
    printf("  Theoretical M           : %.2f\n", (double)M_direct);
    printf("  Theoretical sigma       : %.9f\n", sigma_dir);
    printf("  2-sigma interval        : [%.9f, %.9f]\n",
           rms_direct - two_sigma_dir, rms_direct + two_sigma_dir);
    printf("  3-sigma interval        : [%.9f, %.9f]\n",
           rms_direct - three_sigma_dir, rms_direct + three_sigma_dir);
    printf("  Empirical mean          : %.9f\n", mean_dir);
    printf("  Empirical range         : [%.9f, %.9f]\n", min_dir, max_dir);
    printf("  PASS rate (2-sigma)     : %d / %u (%.1f%%)\n",
           pass_2_dir, repetitions, (100.0 * pass_2_dir) / repetitions);
    printf("  PASS rate (3-sigma)     : %d / %u (%.1f%%)\n",
           pass_3_dir, repetitions, (100.0 * pass_3_dir) / repetitions);

    printf("\nPair mode:\n");
    printf("  Theoretical |cor|       : %.9f\n", aver_rho_sq);
    printf("  Equations needed (M)    : %.2f\n", (double)M_pair);
    printf("  Keystream records       : %llu\n", (unsigned long long)num_symbols_pair);
    printf("  Empirical mean          : %.9f\n", mean_pair);
    printf("  Empirical range         : [%.9f, %.9f]\n", min_pair, max_pair);
    printf("  PASS rate (2-sigma)     : %d / %u (%.1f%%)\n",
           pass_2_pair, repetitions, (100.0 * pass_2_pair) / repetitions);
    printf("  PASS rate (3-sigma)     : %d / %u (%.1f%%)\n",
           pass_3_pair, repetitions, (100.0 * pass_3_pair) / repetitions);

    // Verdict: at least 95% of the repetitions inside the 3-sigma interval
    bool verdict_dir = (pass_3_dir >= (int)(0.95 * repetitions));
    bool verdict_pair = (pass_3_pair >= (int)(0.95 * repetitions));

    printf("\nVerdict (3-sigma): ");
    if (verdict_dir && verdict_pair)
        printf("PASS (both modes within tolerance)\n");
    else if (verdict_dir && !verdict_pair)
        printf("PARTIAL PASS (direct OK, collision slightly marginal)\n");
    else if (!verdict_dir && verdict_pair)
        printf("PARTIAL PASS (collision OK, direct slightly marginal)\n");
    else
        printf("FAIL (significant deviation)\n");

    return 0;
}