/**
 * CFCA initial-state recovery tool with repetition support
 * 
 * This program recovers the 64-bit initial LFSR state of the CFCA stream cipher
 * given a known portion of the state (low bits unknown).
 * 
 * Two modes:
 *   - direct: builds linear equations directly from observed keystream
 *   - pair:   uses collision pairs to reduce unknown bits
 * 
 * Modified: Supports multiple independent runs (--repetitions N) to evaluate
 * success rate and average performance.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <random>

// ----------------------------------------------------------------------
// Bit manipulation helpers with BMI2 support
// ----------------------------------------------------------------------

static inline uint64_t pack_bits(uint64_t value, uint64_t mask) {
#if defined(__BMI2__)
    return __builtin_ia32_pext_di(value, mask);
#else
    uint64_t packed = 0;
    unsigned dst = 0;
    for (unsigned src = 0; src < 64; ++src) {
        if (((mask >> src) & 1U) == 0) continue;
        packed |= ((value >> src) & 1U) << dst;
        ++dst;
    }
    return packed;
#endif
}

static inline uint64_t unpack_bits(uint64_t packed, uint64_t mask) {
#if defined(__BMI2__)
    return __builtin_ia32_pdep_di(packed, mask);
#else
    uint64_t value = 0;
    unsigned src = 0;
    for (unsigned dst = 0; dst < 64; ++dst) {
        if (((mask >> dst) & 1U) == 0) continue;
        value |= ((packed >> src) & 1U) << dst;
        ++src;
    }
    return value;
#endif
}

static const uint64_t DEFAULT_SEED = UINT64_C(0x123456789ABCDEF0);
static const double DEFAULT_SUCCESS_PROBABILITY = 0.999992;
static const unsigned DEFAULT_RECOVER_BITS = 16;
static const unsigned DEFAULT_COLLISION_BITS = 16;
static const unsigned T1_SHIFT = 32;
static const unsigned T2_SHIFT = 0;
static const uint64_t LFSR_REDUCTION = UINT64_C(0x1B);

// ----------------------------------------------------------------------
// GF(2) and GF(16) primitives
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

static inline unsigned count_bits64(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_popcountll(value));
#else
    unsigned count = 0;
    while (value) { value &= value - 1; ++count; }
    return count;
#endif
}

static inline unsigned dot16(uint16_t mask, uint16_t value) {
    return parity64(static_cast<uint64_t>(mask & value));
}

static uint64_t low_bits_mask(unsigned bits) {
    if (bits == 0) return 0;
    if (bits >= 64) return UINT64_MAX;
    return (UINT64_C(1) << bits) - 1;
}

static uint64_t select_high_bits(uint64_t allowed_mask, unsigned number_of_bits) {
    if (count_bits64(allowed_mask) < number_of_bits)
        throw std::runtime_error("not enough positions for the collision key");
    uint64_t result = 0;
    for (int bit = 63; bit >= 0 && number_of_bits; --bit) {
        if (((allowed_mask >> bit) & 1U) != 0) {
            result |= UINT64_C(1) << bit;
            --number_of_bits;
        }
    }
    return result;
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

static uint64_t gf64_multiply(uint64_t a, uint64_t b) {
    uint64_t p = 0;
    while (b) {
        if (b & 1U) p ^= a;
        b >>= 1;
        a = lfsr_clock_bit_value(a);
    }
    return p;
}

static uint64_t gf64_power(uint64_t base, uint64_t exp) {
    uint64_t r = 1;
    while (exp) {
        if (exp & 1U) r = gf64_multiply(r, base);
        exp >>= 1;
        if (exp) base = gf64_multiply(base, base);
    }
    return r;
}

static bool verify_maximal_period() {
    static const uint64_t factors[] = {3,5,17,257,641,65537,6700417};
    const uint64_t order = UINT64_MAX;
    const uint64_t x = 2;
    if (gf64_power(x, order) != 1) return false;
    for (unsigned i = 0; i < sizeof(factors)/sizeof(factors[0]); ++i)
        if (gf64_power(x, order / factors[i]) == 1) return false;
    return true;
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

static uint64_t symbolic_tap_mask(const SymbolicLfsr &s, unsigned sh, uint16_t mask) {
    uint64_t r = 0;
    for (unsigned i = 0; i < 16; ++i)
        if ((mask >> i) & 1U) r ^= s.row[sh + i];
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

static uint64_t apply_symbolic_mask(const std::array<uint64_t, 16> &masks, uint16_t mask) {
    uint64_t r = 0;
    for (unsigned i = 0; i < 16; ++i) if ((mask >> i) & 1U) r ^= masks[i];
    return r;
}

// SplitMix64 PRNG
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
// Automaton8 correlation function (from provided C code)
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
// Adaptive approximation table (pre-computed via FWHT)
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
// Normal CDF inverse
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
// Equation building
// ----------------------------------------------------------------------

struct Equation {
    uint64_t state_mask;
    unsigned observed_bit;
    double corr_abs;
};

static Equation build_equation(const Sample &prev, const Sample &cur, const Sample &next,
                               const std::array<AdaptiveEntry, 16> &table) {
    const uint16_t beta = 0x0008, gamma = 0x0041, m = 0x0008, n = 0x0061, h = 0x0041;
    unsigned u = prev.z & 0x0FU;
    const AdaptiveEntry &e = table[u];
    Equation eq;
    eq.observed_bit = dot16(beta, cur.z) ^ dot16(gamma, next.z) ^ e.tau;
    eq.state_mask = apply_symbolic_mask(prev.T1_bit_masks, e.l) ^
                    apply_symbolic_mask(cur.T1_bit_masks, m) ^
                    apply_symbolic_mask(next.T1_bit_masks, n) ^
                    apply_symbolic_mask(cur.T2_bit_masks, h);
    eq.corr_abs = e.absolute_correlation;
    return eq;
}

// ----------------------------------------------------------------------
// Equation collection
// ----------------------------------------------------------------------

enum RecoveryMode { MODE_DIRECT, MODE_PAIR };

struct PendingEquation {
    uint64_t state_mask;
    unsigned observed_bit;
    double corr_abs;
};

struct CollectionResult {
    std::vector<int64_t> spectrum;
    uint64_t symbols, usable_equations, zero_equations;
    double true_correlation_sum;
};

static CollectionResult collect_equations(
    Cipher &cipher,
    const std::array<AdaptiveEntry, 16> &table,
    RecoveryMode mode,
    uint64_t target_initial_state,
    uint64_t known_state,
    uint64_t coeff_mask,
    uint64_t collision_key_mask,
    unsigned collision_bits,
    uint64_t symbol_limit,
    uint64_t equation_limit)
{
    unsigned coeff_bits = count_bits64(coeff_mask);
    std::size_t cand_count = static_cast<std::size_t>(1) << coeff_bits;
    CollectionResult res;
    res.spectrum.assign(cand_count, 0);
    res.symbols = res.usable_equations = res.zero_equations = 0;
    res.true_correlation_sum = 0.0;

    if (mode == MODE_DIRECT) {
        Sample prev = generate_sample(cipher);
        Sample cur  = generate_sample(cipher);
        Sample next = generate_sample(cipher);

        uint64_t report_interval = symbol_limit / 10;
        if (report_interval == 0) report_interval = 1;

        for (uint64_t idx = 0; idx < symbol_limit && res.usable_equations < equation_limit; ++idx) {
            Equation eq = build_equation(prev, cur, next, table);
            uint64_t coeff = eq.state_mask & coeff_mask;
            uint64_t packed = pack_bits(coeff, coeff_mask);
            unsigned obs = eq.observed_bit ^ parity64(eq.state_mask & known_state);
            res.spectrum[static_cast<std::size_t>(packed)] += obs ? -1 : 1;

            unsigned true_bit = parity64(eq.state_mask & target_initial_state);
            double w = 2.0 * eq.corr_abs;
            res.true_correlation_sum += (eq.observed_bit ^ true_bit) ? -w : w;
            ++res.usable_equations;
            ++res.symbols;

            prev = cur; cur = next; next = generate_sample(cipher);
        }
    } else { // MODE_PAIR
        size_t num_buckets = size_t(1) << collision_bits;
        std::vector<std::vector<PendingEquation>> buckets(num_buckets);

        Sample prev = generate_sample(cipher);
        Sample cur  = generate_sample(cipher);
        Sample next = generate_sample(cipher);

        for (uint64_t idx = 0; idx < symbol_limit; ++idx) {
            Equation eq = build_equation(prev, cur, next, table);
            uint64_t packed_key = pack_bits(eq.state_mask, collision_key_mask);
            PendingEquation peq;
            peq.state_mask = eq.state_mask;
            peq.observed_bit = eq.observed_bit;
            peq.corr_abs = eq.corr_abs;
            buckets[packed_key].push_back(peq);
            ++res.symbols;

            prev = cur; cur = next; next = generate_sample(cipher);
        }

        for (auto& bucket : buckets) {
            size_t m = bucket.size();
            for (size_t i = 0; i < m && res.usable_equations < equation_limit; ++i) {
                for (size_t j = i + 1; j < m && res.usable_equations < equation_limit; ++j) {
                    const PendingEquation& first = bucket[i];
                    const PendingEquation& second = bucket[j];

                    uint64_t delta_mask = first.state_mask ^ second.state_mask;
                    unsigned pair_obs    = first.observed_bit ^ second.observed_bit;

                    if ((delta_mask & collision_key_mask) != 0)
                        throw std::logic_error("collision key did not cancel");

                    uint64_t coeff = delta_mask & coeff_mask;
                    if (coeff == 0) {
                        ++res.zero_equations;
                    } else {
                        unsigned obs = pair_obs ^ parity64(delta_mask & known_state);
                        uint64_t packed = pack_bits(coeff, coeff_mask);
                        res.spectrum[static_cast<std::size_t>(packed)] += obs ? -1 : 1;

                        unsigned true_bit = parity64(delta_mask & target_initial_state);
                        double w = 2.0 * (first.corr_abs * second.corr_abs);
                        res.true_correlation_sum += (pair_obs ^ true_bit) ? -w : w;
                        ++res.usable_equations;
                    }
                }
            }
            if (res.usable_equations >= equation_limit) break;
        }
    }
    return res;
}

// ----------------------------------------------------------------------
// FWHT and recovery
// ----------------------------------------------------------------------

static void fwht(std::vector<int64_t> &vals) {
    std::size_t size = vals.size();
    for (std::size_t len = 1; len < size; len <<= 1)
        for (std::size_t st = 0; st < size; st += len << 1)
            for (std::size_t o = 0; o < len; ++o) {
                int64_t a = vals[st + o], b = vals[st + o + len];
                vals[st + o] = a + b;
                vals[st + o + len] = a - b;
            }
}

struct RecoveryResult {
    uint64_t recovered_state;
    uint64_t recovered_unknown_packed;
    int64_t best_score, true_score;
    std::size_t true_rank;
    bool success;
};

static RecoveryResult recover_state_direct(
    std::vector<int64_t> spectrum,
    uint64_t coeff_mask,
    uint64_t known_state,
    uint64_t true_state_full)
{
    fwht(spectrum);
    std::size_t best_idx = 0;
    int64_t best = spectrum[0];
    for (std::size_t i = 1; i < spectrum.size(); ++i)
        if (spectrum[i] > best) { best = spectrum[i]; best_idx = i; }

    uint64_t true_idx = pack_bits(true_state_full, coeff_mask);
    int64_t true_sc = spectrum[static_cast<std::size_t>(true_idx)];
    std::size_t rank = 1;
    for (std::size_t i = 0; i < spectrum.size(); ++i)
        if (spectrum[i] > true_sc) ++rank;

    uint64_t recovered_unknown = unpack_bits(static_cast<uint64_t>(best_idx), coeff_mask);
    RecoveryResult r;
    r.recovered_state = known_state | recovered_unknown;
    r.recovered_unknown_packed = static_cast<uint64_t>(best_idx);
    r.best_score = best;
    r.true_score = true_sc;
    r.true_rank = rank;
    r.success = (r.recovered_state == true_state_full);
    return r;
}

// ----------------------------------------------------------------------
// Helper for pair mode symbol count estimation
// ----------------------------------------------------------------------

struct VerifySample { uint16_t T2; uint16_t z; };

static long double expected_disjoint_pairs(long double records, long double buckets) {
    if (records <= 0.0L) return 0.0L;
    if (buckets <= 1.0L) return std::floor(records / 2.0L);
    long double logt = std::log1pl(-2.0L / buckets);
    long double odd  = 1.0L - std::exp(records * logt);
    return records / 2.0L - buckets * odd / 4.0L;
}

static long double required_collision_records(long double target_pairs, unsigned col_bits) {
    long double buckets = std::ldexp(1.0L, static_cast<int>(col_bits));
    long double low = 0.0L, high = std::max(2.0L * target_pairs + buckets, 1.0L);
    while (expected_disjoint_pairs(high, buckets) < target_pairs) high *= 2.0L;
    for (unsigned i = 0; i < 100; ++i) {
        long double mid = (low + high) / 2.0L;
        if (expected_disjoint_pairs(mid, buckets) < target_pairs) low = mid;
        else high = mid;
    }
    return high;
}

// ----------------------------------------------------------------------
// Options and main
// ----------------------------------------------------------------------

struct Options {
    RecoveryMode mode;
    unsigned recover_bits, collision_bits;
    uint64_t symbol_limit, equation_limit, seed;
    double success_probability;
    bool seed_provided;
    unsigned repetitions;
};

static void print_usage(const char *prog) {
    std::printf("Usage: %s [options]\n"
                "  --mode direct|pair  Recovery mode (default direct)\n"
                "  --B N               Unknown low state bits (default 16)\n"
                "  --K N               Pair collision-key bits (default 16, must be < B)\n"
                "  --symbols N         Max symbols to process (default: theoretical M for direct, sqrt(M*2^(K+1)) for pair)\n"
                "  --equations N       Max usable equations (default: theoretical M)\n"
                "  --ps P              Target success probability (default 0.999992)\n"
                "  --seed N            Experiment seed (if not given, random)\n"
                "  --repetitions N     Number of independent runs (default 100)\n"
                "  --help              Show this message\n", prog);
}

static Options parse_options(int argc, char **argv) {
    Options opt;
    opt.mode = MODE_DIRECT;
    opt.recover_bits = DEFAULT_RECOVER_BITS;
    opt.collision_bits = DEFAULT_COLLISION_BITS;
    opt.symbol_limit = 0;   // 0 means auto
    opt.equation_limit = 0;  // 0 means auto
    opt.seed = DEFAULT_SEED;
    opt.success_probability = DEFAULT_SUCCESS_PROBABILITY;
    opt.seed_provided = false;
    opt.repetitions = 1;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--mode") == 0 && i+1 < argc) {
            const char *v = argv[++i];
            if (std::strcmp(v, "direct") == 0) opt.mode = MODE_DIRECT;
            else if (std::strcmp(v, "pair") == 0) opt.mode = MODE_PAIR;
            else throw std::runtime_error("--mode must be direct or pair");
        } else if (std::strcmp(argv[i], "--B") == 0 && i+1 < argc)
            opt.recover_bits = static_cast<unsigned>(std::strtoul(argv[++i], NULL, 0));
        else if (std::strcmp(argv[i], "--K") == 0 && i+1 < argc)
            opt.collision_bits = static_cast<unsigned>(std::strtoul(argv[++i], NULL, 0));
        else if (std::strcmp(argv[i], "--symbols") == 0 && i+1 < argc)
            opt.symbol_limit = std::strtoull(argv[++i], NULL, 0);
        else if (std::strcmp(argv[i], "--equations") == 0 && i+1 < argc)
            opt.equation_limit = std::strtoull(argv[++i], NULL, 0);
        else if (std::strcmp(argv[i], "--ps") == 0 && i+1 < argc)
            opt.success_probability = std::strtod(argv[++i], NULL);
        else if (std::strcmp(argv[i], "--seed") == 0 && i+1 < argc) {
            opt.seed = std::strtoull(argv[++i], NULL, 0);
            opt.seed_provided = true;
        } else if (std::strcmp(argv[i], "--repetitions") == 0 && i+1 < argc) {
            opt.repetitions = static_cast<unsigned>(std::strtoul(argv[++i], NULL, 0));
            if (opt.repetitions == 0) opt.repetitions = 1;
        } else if (std::strcmp(argv[i], "--help") == 0) { print_usage(argv[0]); std::exit(0); }
        else throw std::runtime_error("invalid option");
    }
    return opt;
}

// ----------------------------------------------------------------------
// Main (with repetition support)
// ----------------------------------------------------------------------

int main(int argc, char **argv) {
    try {
        Options opt = parse_options(argc, argv);

        if (opt.recover_bits == 0 || opt.recover_bits > 24)
            throw std::runtime_error("B must be between 1 and 24");
        if (opt.collision_bits == 0 || opt.collision_bits > 24)
            throw std::runtime_error("K must be between 1 and 24");
        if (opt.mode == MODE_PAIR && opt.collision_bits >= opt.recover_bits)
            throw std::runtime_error("K must be strictly less than B in pair mode");
        if (!(opt.success_probability > 0.0 && opt.success_probability < 1.0))
            throw std::runtime_error("success probability must be in (0,1)");

        // ---- Compute adaptive table (same for all repetitions) ----
        const uint16_t beta_fixed = 0x0008;
        const uint16_t gamma_fixed = 0x0041;
        const uint16_t h_fixed = 0x0041;
        const uint16_t n_fixed = 0x0061;

        std::vector<int> xi_list = {
            0x17, 0x28, 0x3F, 0x41, 0x56, 0x69, 0x7E,
            0x83, 0x94, 0xAB, 0xBC, 0xC2, 0xD5, 0xEA, 0xFD
        };

        auto table = compute_adaptive_table(xi_list, beta_fixed, gamma_fixed, h_fixed, n_fixed);

        double aver_rho_sq = 0.0;
        for (auto &e : table) aver_rho_sq += e.absolute_correlation * e.absolute_correlation;
        aver_rho_sq /= 16.0;
        double rms_direct = std::sqrt(aver_rho_sq);

        // ---- Theoretical M ----
        double decoder_rho_sq;
        if (opt.mode == MODE_DIRECT)
            decoder_rho_sq = aver_rho_sq;
        else
            decoder_rho_sq = aver_rho_sq * aver_rho_sq;

        unsigned effective_bits = (opt.mode == MODE_DIRECT) ? opt.recover_bits : opt.recover_bits - opt.collision_bits;
        double z = inverse_normal_cdf(opt.success_probability);
        // Selcuk formula with B (not B+1) in the sqrt term: sqrt(2 * B * ln2)
        double sqrt_term = std::sqrt(2.0 * std::log(2.0) * effective_bits);
        long double formula_M = (1.0L / (long double)decoder_rho_sq) *
                                std::pow((long double)(z + sqrt_term), 2);
        if (!std::isfinite(static_cast<double>(formula_M)) || formula_M > static_cast<long double>(UINT64_MAX))
            throw std::runtime_error("calculated M does not fit in uint64_t");
        uint64_t eq_limit = static_cast<uint64_t>(std::ceil(formula_M));
        
        // ---- Determine symbol limit (AUTO) ----
        uint64_t sym_limit = opt.symbol_limit;
        if (sym_limit == 0) {
            if (opt.mode == MODE_DIRECT) {
                sym_limit = eq_limit;
            } else {
                // N = sqrt(M * 2^(K+1))
                long double need_sym = std::sqrt((long double)formula_M *
                                                 std::ldexp(1.0L, opt.collision_bits + 1));
                sym_limit = static_cast<uint64_t>(std::ceil(need_sym));
            }
        }
        
        // Use theoretical M as equation limit if not set
        uint64_t eq_limit_final = (opt.equation_limit == 0) ? eq_limit : opt.equation_limit;

        std::printf("===== CFCA initial-state recovery (with repetitions) =====\n");
        std::printf("Mode: %s\n", opt.mode == MODE_DIRECT ? "direct" : "pair");
        std::printf("B = %u, K = %u\n", opt.recover_bits, opt.collision_bits);
        std::printf("Target success prob = %.6f\n", opt.success_probability);
        std::printf("Theoretical M = %llu\n", (unsigned long long)eq_limit);
        std::printf("Symbol limit = %llu\n", (unsigned long long)sym_limit);
        std::printf("Equation limit = %llu\n", (unsigned long long)eq_limit_final);
        std::printf("Repetitions = %u\n", opt.repetitions);
        std::printf("aver_rho_sq = %.9f\n", aver_rho_sq);
        std::printf("========================================================\n");

        // ---- Loop over repetitions ----
        unsigned successes = 0;
        std::vector<uint64_t> eq_used_list;
        std::vector<std::size_t> rank_list;

        for (unsigned rep = 0; rep < opt.repetitions; ++rep) {
            uint64_t current_seed;
            if (opt.seed_provided) {
                current_seed = opt.seed + rep;
            } else {
                std::random_device rd;
                current_seed = (static_cast<uint64_t>(rd()) << 32) | rd();
                if (current_seed == 0) current_seed = 1;
            }

            Cipher cipher = initialize_cipher(current_seed);
            uint64_t target_initial = get_state(cipher.lfsr);
            uint64_t unknown_mask = low_bits_mask(opt.recover_bits);
            uint64_t known_state = target_initial & ~unknown_mask;

            uint64_t collision_key_mask = 0;
            uint64_t coeff_mask = unknown_mask;
            if (opt.mode == MODE_PAIR) {
                collision_key_mask = select_high_bits(unknown_mask, opt.collision_bits);
                coeff_mask = unknown_mask & ~collision_key_mask;
            }

            // Capture verification samples
            std::vector<VerifySample> verify_samples;
            uint16_t init_R1 = cipher.R1, init_R2 = cipher.R2, init_R3 = cipher.R3;
            for (int i = 0; i < 5; ++i) {
                Sample s = generate_sample(cipher);
                verify_samples.push_back({s.T2, s.z});
            }

            CollectionResult collection = collect_equations(
                cipher, table, opt.mode, target_initial, known_state, coeff_mask,
                collision_key_mask, opt.collision_bits, sym_limit, eq_limit_final);

            if (collection.usable_equations == 0) {
                std::printf("Rep %u: seed=%016llX, no equations collected.\n",
                            rep, (unsigned long long)current_seed);
                continue;
            }

            RecoveryResult stage1 = recover_state_direct(collection.spectrum, coeff_mask,
                                                         known_state, target_initial);

            bool final_success = false;
            uint64_t final_recovered_state = 0;

            if (opt.mode == MODE_DIRECT) {
                final_recovered_state = stage1.recovered_state;
                final_success = stage1.success;
            } else { // MODE_PAIR
                uint64_t partial_unknown = stage1.recovered_state ^ known_state;
                uint64_t max_guess = UINT64_C(1) << opt.collision_bits;
                bool found = false;
                uint64_t found_state = 0;

                for (uint64_t guess = 0; guess < max_guess; ++guess) {
                    uint64_t collision_part = unpack_bits(guess, collision_key_mask);
                    uint64_t candidate = known_state | partial_unknown | collision_part;

                    Cipher test_cipher;
                    set_state(test_cipher.lfsr, candidate);
                    test_cipher.R1 = init_R1;
                    test_cipher.R2 = init_R2;
                    test_cipher.R3 = init_R3;

                    bool match = true;
                    for (const auto &v : verify_samples) {
                        uint16_t T1 = get_T1(test_cipher.lfsr);
                        uint16_t T2 = get_T2(test_cipher.lfsr);
                        uint16_t z_out = static_cast<uint16_t>(static_cast<uint16_t>(test_cipher.R1 + T1) ^ test_cipher.R2);
                        if (T2 != v.T2 || z_out != v.z) { match = false; break; }
                        fsm_update(test_cipher, T2);
                        lfsr_clock_word(test_cipher.lfsr);
                    }
                    if (match) {
                        found = true;
                        found_state = candidate;
                        break;
                    }
                }
                if (found) {
                    final_recovered_state = found_state;
                    final_success = true;
                }
            }

            if (final_success) {
                ++successes;
                eq_used_list.push_back(collection.usable_equations);
                rank_list.push_back(stage1.true_rank);
            }

            std::printf("Rep %u: seed=%016llX, eq=%llu, rank=%llu/%llu, %s\n",
                        rep, (unsigned long long)current_seed,
                        (unsigned long long)collection.usable_equations,
                        (unsigned long long)stage1.true_rank,
                        (unsigned long long)(UINT64_C(1) << effective_bits),
                        final_success ? "SUCCESS" : "FAIL");
        }

        // ---- Summary ----
        double success_rate = (opt.repetitions > 0) ? (double)successes / opt.repetitions : 0.0;
        std::printf("\n===== Summary =====\n");
        std::printf("Total runs: %u\n", opt.repetitions);
        std::printf("Successes: %u\n", successes);
        std::printf("Success rate: %.1f%%\n", success_rate * 100.0);

        if (!eq_used_list.empty()) {
            double avg_eq = 0.0;
            for (auto v : eq_used_list) avg_eq += (double)v;
            avg_eq /= eq_used_list.size();
            double avg_rank = 0.0;
            for (auto v : rank_list) avg_rank += (double)v;
            avg_rank /= rank_list.size();
            std::printf("Average equations used (successful): %.0f\n", avg_eq);
            std::printf("Average rank (successful): %.1f\n", avg_rank);
        }

        // Compare to target
        double target_rate = opt.success_probability;
        bool success_threshold_met = (success_rate >= target_rate);
        if (success_threshold_met) {
            std::printf("\nSUCCESS: Success rate %.1f%% meets target %.1f%%\n", 
                        success_rate * 100.0, target_rate * 100.0);
        } else {
            std::printf("\nFAIL: Success rate %.1f%% is below target %.1f%%\n", 
                        success_rate * 100.0, target_rate * 100.0);
        }

        return success_threshold_met ? EXIT_SUCCESS : 2; 
    } catch (const std::bad_alloc &) {
        std::fprintf(stderr, "Error: insufficient memory. Reduce B or K.\n");
        return EXIT_FAILURE;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return EXIT_FAILURE;
    }
}