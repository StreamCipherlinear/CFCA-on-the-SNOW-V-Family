/**
 * verify_conditional_independence_KL.cpp
 *
 * This program verifies whether the conditional approximation of e2
 * (under condition u = (invR2 + T1_prev) ^ invR3) is independent of e3.
 *
 * It builds an adaptive linear approximation for e2 in the form:
 *   e2_cond = parity(l * T1_prev) ^ tau ^ parity(xi * R2) ^ parity(h * R3)
 * where (l, tau) are chosen by Walsh-Hadamard transform over a merged set
 * of samples for all given xi values.
 *
 * Then, for each xi in the list, it computes the KL divergence between
 * the joint distribution of (e2_cond, e3) and the product of marginals.
 *
 * Usage: ./verify_KL (no arguments, parameters fixed)
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdlib>

// ============================================================
// GF(2^4) and S-box definitions (same as Q1)
// ============================================================

static uint8_t gf4_mul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    while (b) {
        if (b & 1) p ^= a;
        a <<= 1;
        if (a & 0x10) a ^= 0x13;
        b >>= 1;
    }
    return p & 0xF;
}

static const uint8_t SBOX[16] = {
    0x9, 0x4, 0xA, 0xB, 0xD, 0x1, 0x8, 0x5,
    0x6, 0x2, 0x0, 0x3, 0xC, 0xE, 0xF, 0x7
};

static inline uint8_t sr(uint8_t x) {
    return SBOX[x & 0xF];
}

static uint8_t S8_transform(uint8_t x) {
    uint8_t x0 = x & 0xF;
    uint8_t x1 = (x >> 4) & 0xF;
    uint8_t y0 = sr(x0);
    uint8_t y1 = sr(x1);
    uint8_t z0 = y0 ^ gf4_mul(4, y1);
    uint8_t z1 = gf4_mul(4, y0) ^ y1;
    return (z0 & 0xF) | ((z1 & 0xF) << 4);
}

// ============================================================
// Original e2 and e3 formulas (same as Q1)
// ============================================================

static inline uint8_t compute_e2(uint8_t invR2, uint8_t invR3, uint8_t R2, uint8_t R3,
                                 uint8_t T1_prev, uint8_t alpha4, uint8_t l8,
                                 uint8_t xi8, uint8_t h8) {
    uint8_t add = (invR2 + T1_prev) & 0xFF;
    int p1 = __builtin_parity((unsigned int)(alpha4 & add));
    int p2 = __builtin_parity((unsigned int)(l8 & T1_prev));
    int p3 = __builtin_parity((unsigned int)(xi8 & R2));
    int p4 = __builtin_parity((unsigned int)(alpha4 & invR3));
    int p5 = __builtin_parity((unsigned int)(h8 & R3));
    return (p1 ^ p2 ^ p3 ^ p4 ^ p5) & 0x1;
}

static inline uint8_t compute_e3(uint8_t R2, uint8_t R3, uint8_t T2, uint8_t T1_next,
                                 uint8_t gamma8, uint8_t n8, uint8_t beta4,
                                 uint8_t xi8, uint8_t h8) {
    uint8_t inner = R3 ^ T2;
    uint8_t xored = (inner + R2) & 0xFF;
    uint8_t s_out = (xored + T1_next) & 0xFF;
    int p1 = __builtin_parity((unsigned int)(gamma8 & s_out));
    int p2 = __builtin_parity((unsigned int)(n8 & T1_next));
    int p3 = __builtin_parity((unsigned int)((beta4 ^ xi8) & R2));
    int p4 = __builtin_parity((unsigned int)(h8 & (R3 ^ T2)));
    return (p1 ^ p2 ^ p3 ^ p4) & 0x1;
}

// ============================================================
// 8-bit Walsh-Hadamard Transform
// ============================================================

static void fwht_8bit(std::vector<int>& data) {
    int n = (int)data.size();  // must be 256
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; ++j) {
                int a = data[i + j];
                int b = data[i + j + len];
                data[i + j] = a + b;
                data[i + j + len] = a - b;
            }
        }
    }
}

// ============================================================
// Adaptive table entry: (l, tau) for each condition u
// ============================================================

struct AdaptiveEntry {
    uint8_t l;      // optimal linear mask (8-bit)
    uint8_t tau;    // bias (0 or 1)
};

// ============================================================
// Build adaptive table by merging samples over all given xi values
// Condition: u = (invR2 + T1_prev) ^ invR3
// For each u, we accumulate the sign of e2 over all (invR2, invR3, T1_prev)
// and all xi in the list, then perform 8-bit WHT to select the best l and tau.
// ============================================================

static std::array<AdaptiveEntry, 256> build_adaptive_table_merged(
        const std::vector<uint8_t>& xi_list,
        uint8_t alpha4,
        uint8_t l8,
        uint8_t h8) {

    std::array<AdaptiveEntry, 256> table;
    const int STATE_SIZE = 1 << 8;  // 256 possible values of T1_prev

    // spectra[u][T1_prev] accumulates +1 if e2=0, -1 if e2=1
    std::vector<std::vector<int>> spectra(256, std::vector<int>(STATE_SIZE, 0));

    #pragma omp parallel for schedule(dynamic)
    for (uint32_t rr = 0; rr < 65536; ++rr) {
        uint8_t invR2 = rr & 0xFF;
        uint8_t invR3 = (rr >> 8) & 0xFF;
        uint8_t R2 = S8_transform(invR2);
        uint8_t R3 = S8_transform(invR3);

        for (uint32_t t = 0; t < 256; ++t) {
            uint8_t T1_prev = (uint8_t)t;
            uint8_t add = (invR2 + T1_prev) & 0xFF;
            uint8_t u = add ^ invR3;

            // For each xi in the list, compute e2 and add its sign
            for (uint8_t xi : xi_list) {
                uint8_t e2 = compute_e2(invR2, invR3, R2, R3, T1_prev,
                                        alpha4, l8, xi, h8);
                int val = (e2 == 0) ? 1 : -1;
                #pragma omp atomic
                spectra[u][T1_prev] += val;
            }
        }
    }

    // For each u, run WHT and pick the mask with largest absolute correlation
    for (int u = 0; u < 256; ++u) {
        std::vector<int>& spec = spectra[u];
        fwht_8bit(spec);  // after transform, spec[l] = sum_T sgn(T) * (-1)^{dot(l,T)}

        int best_l = 0;
        int best_val = spec[0];
        for (int l = 1; l < STATE_SIZE; ++l) {
            if (std::abs(spec[l]) > std::abs(best_val)) {
                best_val = spec[l];
                best_l = l;
            }
        }
        table[u].l = (uint8_t)best_l;
        // if best_val is negative, we flip the output with tau=1
        table[u].tau = (best_val < 0) ? 1 : 0;
    }

    return table;
}

// ============================================================
// Evaluate KL divergence for a single xi using a fixed table
// The approximated e2 is:
//   e2_cond = parity(l * T1_prev) ^ tau ^ parity(xi * R2) ^ parity(h * R3)
// where (l, tau) are taken from the table based on u.
// The joint distribution of (e2_cond, e3) is computed exactly by
// enumerating all inputs (invR2, invR3, T1_prev, T1_next, T2).
// ============================================================

static double evaluate_KL_with_fixed_table(
        const std::array<AdaptiveEntry, 256>& table,
        uint8_t xi8,
        uint8_t alpha4,
        uint8_t beta4,
        uint8_t gamma8,
        uint8_t l8,
        uint8_t n8,
        uint8_t h8) {

    uint64_t j00 = 0, j01 = 0, j10 = 0, j11 = 0;

    #pragma omp parallel for reduction(+:j00, j01, j10, j11) schedule(dynamic)
    for (uint32_t rr = 0; rr < 65536; ++rr) {
        uint8_t invR2 = rr & 0xFF;
        uint8_t invR3 = (rr >> 8) & 0xFF;
        uint8_t R2 = S8_transform(invR2);
        uint8_t R3 = S8_transform(invR3);

        // Count frequencies of e2_cond (conditional part) over T1_prev
        uint64_t cnt2[2] = {0, 0};
        for (uint32_t t = 0; t < 256; ++t) {
            uint8_t T1_prev = (uint8_t)t;
            uint8_t add = (invR2 + T1_prev) & 0xFF;
            uint8_t u = add ^ invR3;
            const AdaptiveEntry& e = table[u];

            int parity_l = __builtin_parity((unsigned int)(e.l & T1_prev));
            int p3 = __builtin_parity((unsigned int)(xi8 & R2));
            int p5 = __builtin_parity((unsigned int)(h8 & R3));
            uint8_t n2_cond = (parity_l ^ e.tau ^ p3 ^ p5) & 0x1;
            cnt2[n2_cond]++;
        }

        // Count frequencies of e3 over (T1_next, T2)
        uint64_t cnt3[2] = {0, 0};
        for (uint32_t tt = 0; tt < 65536; ++tt) {
            uint8_t T1_next = (uint8_t)(tt & 0xFF);
            uint8_t T2 = (uint8_t)((tt >> 8) & 0xFF);
            uint8_t e3 = compute_e3(R2, R3, T2, T1_next, gamma8, n8, beta4, xi8, h8);
            cnt3[e3]++;
        }

        // Accumulate joint counts (conditional independence inside the loop)
        j00 += cnt2[0] * cnt3[0];
        j01 += cnt2[0] * cnt3[1];
        j10 += cnt2[1] * cnt3[0];
        j11 += cnt2[1] * cnt3[1];
    }

    uint64_t joint[2][2] = {{j00, j01}, {j10, j11}};
    double total = (double)(j00 + j01 + j10 + j11);
    if (total == 0.0) return 0.0;

    // Empirical joint distribution
    double P[2][2];
    double marg_e2[2] = {0.0, 0.0};
    double marg_e3[2] = {0.0, 0.0};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            P[i][j] = (double)joint[i][j] / total;
            marg_e2[i] += P[i][j];
            marg_e3[j] += P[i][j];
        }
    }

    // KL divergence D_KL( P(e2_cond, e3) || P(e2_cond)*P(e3) )
    double D_KL = 0.0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            double indep = marg_e2[i] * marg_e3[j];
            if (P[i][j] > 0.0 && indep > 0.0)
                D_KL += P[i][j] * log2(P[i][j] / indep);
        }
    }

    return D_KL;
}

// ============================================================
// Main program
// ============================================================

int main() {
    // Fixed parameters (as in Q1)
    const uint8_t alpha4 = 0xC;
    const uint8_t beta4  = 0x8;
    const uint8_t gamma8 = 0x41;
    const uint8_t l8     = 0x8;   // fixed mask in original e2
    const uint8_t n8     = 0x61;
    const uint8_t h8     = 0x41;

    // List of xi values (as in Q3)
    std::vector<uint8_t> xi_list = {
        0x17, 0x28, 0x3F, 0x41, 0x56, 0x69, 0x7E,
        0x83, 0x94, 0xAB, 0xBC, 0xC2, 0xD5, 0xEA, 0xFD
    };

    printf("Building adaptive table by merging over all xi in the list...\n");
    auto table = build_adaptive_table_merged(xi_list, alpha4, l8, h8);

    printf("\nComputing KL divergence for each xi using the merged table...\n");
    std::vector<double> kl_values;
    kl_values.reserve(xi_list.size());

    double sum_kl = 0.0;
    double max_kl = -1.0, min_kl = 1e100;
    uint8_t max_xi = 0, min_xi = 0;

    for (uint8_t xi : xi_list) {
        double dkl = evaluate_KL_with_fixed_table(table, xi,
                                                  alpha4, beta4, gamma8,
                                                  l8, n8, h8);
        printf("xi = 0x%02X  D_KL = %.8f\n", xi, dkl);
        kl_values.push_back(dkl);
        sum_kl += dkl;
        if (dkl > max_kl) { max_kl = dkl; max_xi = xi; }
        if (dkl < min_kl) { min_kl = dkl; min_xi = xi; }
    }

    double avg_kl = sum_kl / (double)xi_list.size();

    printf("\n========== SUMMARY ==========\n");
    printf("Number of xi values: %zu\n", xi_list.size());
    printf("Average KL divergence: %.8f\n", avg_kl);
    printf("Minimum KL: %.8f  (xi=0x%02X)\n", min_kl, min_xi);
    printf("Maximum KL: %.8f  (xi=0x%02X)\n", max_kl, max_xi);

    return 0;
}