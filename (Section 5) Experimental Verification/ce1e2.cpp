#include <stdio.h>
#include <stdint.h>
#include <math.h>    
#include <stdlib.h>  

// ========== Dot product  ==========
static inline int dot(uint8_t a, uint8_t b) {
    return __builtin_parity((unsigned int)(a & b));
}

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

static inline uint8_t sr(uint8_t x) { return SBOX[x & 0xF]; }

static uint8_t S8_transform(uint8_t x) {
    uint8_t x0 = x & 0xF, x1 = (x >> 4) & 0xF;
    uint8_t y0 = sr(x0), y1 = sr(x1);
    uint8_t z0 = y0 ^ gf4_mul(4, y1);
    uint8_t z1 = gf4_mul(4, y0) ^ y1;
    return (z0 & 0xF) | ((z1 & 0xF) << 4);
}

// Compute the linear correlation (bias) of e1
// Formula: e1 = beta*x ^ m*S8(x) ^ gamma*x
double compute_e1_correlation(uint8_t beta, uint8_t m, uint8_t gamma) {
    int cnt = 0;
    for (int x = 0; x < 16; x++) {
        for (int y = 0; y < 16; y++) {
            int z = (x + y) & 0xf;
            if ((dot(beta, z) ^ dot(m, x) ^ dot(gamma, S8_transform(y))) == 0)
                cnt++;
            else
                cnt--;
        }
    }
    return cnt / 256.0;   
}

int main() {
    // ========== Search for the optimal m that maximizes |cor(e1)|, with beta=0x8, gamma=0x41 ==========
    uint8_t fixed_beta = 0x8;
    uint8_t fixed_gamma = 0x41;   
    uint8_t best_m = 0;
    double best_abs_cor_e1 = 0.0;
    double best_cor_e1 = 0.0;          

    printf("===== Searching best m for e1 (beta=0x%X, gamma=0x%X) =====\n",
           fixed_beta, fixed_gamma);
    for (int m = 0; m < 16; m++) {
        double cor = compute_e1_correlation(fixed_beta, (uint8_t)m, fixed_gamma);
        double abs_cor = fabs(cor);
        if (abs_cor > best_abs_cor_e1) {
            best_abs_cor_e1 = abs_cor;
            best_cor_e1 = cor;         
            best_m = m;
        }
    }
    
    printf("Optimal m = 0x%02X, |cor(e1)| = %.6f, cor(e1) = %+.6f\n\n",
           best_m, best_abs_cor_e1, best_cor_e1);

    // ========== Search for the optimal (alpha, l) for e2 ==========
    uint8_t h  = 0x41;   
    uint8_t xi = 0x69;   

    int best_alpha = 0, best_l = 0;
    double best_abs_cor = 0.0;
    double best_cor = 0.0;            // save the actual total correlation

    printf("Search for best (alpha, l) maximizing |cor(e2)|\n");
    printf("Fixed masks: h=0x%02X, xi=0x%02X\n\n", h, xi);

    for (int alpha = 0; alpha < 16; alpha++) {
        // ---------- Part 1: depends only on invR3 ----------
        int s1 = 0;
        for (int x = 0; x < 16; x++) {
            if ((dot(alpha, x) ^ dot(h, S8_transform(x))) == 0)
                s1++;
            else
                s1--;
        }
        double cor1 = s1 / 16.0;

        // ---------- Part 2: depends on invR2 and T1_prev ----------
        for (int l = 0; l < 16; l++) {
            int s2 = 0;
            for (int x = 0; x < 16; x++) {
                uint8_t R2 = S8_transform(x);
                for (int y = 0; y < 16; y++) {
                    uint8_t z = (x + y) & 0xF;
                    if ((dot(alpha, z) ^ dot(xi, R2) ^ dot(l, y)) == 0)
                        s2++;
                    else
                        s2--;
                }
            }
            double cor2 = s2 / 256.0;
            double cor_total = cor1 * cor2;
            double abs_cor = fabs(cor_total);

            if (abs_cor > best_abs_cor) {
                best_abs_cor = abs_cor;
                best_cor = cor_total;       
                best_alpha = alpha;
                best_l = l;
            }
        }
    }

    printf("Optimal combination found:\n");
    printf("  alpha = 0x%X, l = 0x%X\n", best_alpha, best_l);
    printf("  |cor(e2)| = %.6f, cor(e2) = %+.6f\n", best_abs_cor, best_cor);

    return 0;
}