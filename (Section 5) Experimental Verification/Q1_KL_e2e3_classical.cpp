#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// ==================== GF(2^4) and Small AES S-box ====================
static uint8_t gf4_mul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    while (b) { if (b & 1) p ^= a; a <<= 1; if (a & 0x10) a ^= 0x13; b >>= 1; }
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

// ==================== e2 and e3 formulas ====================
static inline uint8_t compute_e2(uint8_t invR2, uint8_t invR3, uint8_t R2, uint8_t R3, uint8_t T1_prev,
                                 uint8_t alpha4, uint8_t l8, uint8_t xi8, uint8_t h8) {
    uint8_t add = (invR2 + T1_prev) & 0xFF;
    int p1 = __builtin_parity((unsigned int)(alpha4 & add));
    int p2 = __builtin_parity((unsigned int)(l8 & T1_prev));
    int p3 = __builtin_parity((unsigned int)(xi8 & R2));
    int p4 = __builtin_parity((unsigned int)(alpha4 & invR3));
    int p5 = __builtin_parity((unsigned int)(h8 & R3));
    return (p1 ^ p2 ^ p3 ^ p4 ^ p5) & 0x1;
}

static inline uint8_t compute_e3(uint8_t R2, uint8_t R3, uint8_t T2, uint8_t T1_next,
                                 uint8_t gamma8, uint8_t n8, uint8_t beta4, uint8_t xi8, uint8_t h8) {
    uint8_t inner = R3 ^ T2;
    uint8_t xored = (inner + R2)  & 0xFF;
    uint8_t s_out = (xored + T1_next) & 0xFF;
    int p1 = __builtin_parity((unsigned int)(gamma8 & s_out));
    int p2 = __builtin_parity((unsigned int)(n8 & T1_next));
    int p3 = __builtin_parity((unsigned int)((beta4 ^ xi8) & R2));
    int p4 = __builtin_parity((unsigned int)(h8 & (R3 ^ T2)));
    return (p1 ^ p2 ^ p3 ^ p4) & 0x1;
}

// ==================== Evaluation function (returns KL divergence) ====================
double evaluate_independence_like_paper(uint8_t alpha4, uint8_t beta4, uint8_t gamma8,
                                        uint8_t xi8, uint8_t l8, uint8_t n8, uint8_t h8) {
    uint64_t j00 = 0, j01 = 0, j10 = 0, j11 = 0;

#pragma omp parallel for reduction(+:j00, j01, j10, j11) schedule(dynamic)
    for (uint32_t rr = 0; rr < 65536; rr++) {
        uint8_t invR2 = rr & 0xFF;
        uint8_t invR3 = (rr >> 8) & 0xFF;
        uint8_t R2 = S8_transform(invR2);
        uint8_t R3 = S8_transform(invR3);

        uint64_t cnt2[2] = {0, 0};
        for (uint32_t t = 0; t < 256; t++) {
            uint8_t T1_prev = (uint8_t)t;
            uint8_t e2 = compute_e2(invR2, invR3, R2, R3, T1_prev, alpha4, l8, xi8, h8);
            cnt2[e2]++;
        }

        uint64_t cnt3[2] = {0, 0};
        for (uint32_t tt = 0; tt < 65536; tt++) {
            uint8_t T1_next = (uint8_t)(tt & 0xFF);
            uint8_t T2 = (uint8_t)((tt >> 8) & 0xFF);
            uint8_t e3 = compute_e3(R2, R3, T2, T1_next, gamma8, n8, beta4, xi8, h8);
            cnt3[e3]++;
        }

        j00 += cnt2[0] * cnt3[0];
        j01 += cnt2[0] * cnt3[1];
        j10 += cnt2[1] * cnt3[0];
        j11 += cnt2[1] * cnt3[1];
    }

    uint64_t joint[2][2] = {{j00, j01}, {j10, j11}};

    double total = (double)(j00 + j01 + j10 + j11);
    double P[2][2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            P[i][j] = (double)joint[i][j] / total;

    double marg_e2[2] = {P[0][0] + P[0][1], P[1][0] + P[1][1]};
    double marg_e3[2] = {P[0][0] + P[1][0], P[0][1] + P[1][1]};

    double D_KL = 0.0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            double indep = marg_e2[i] * marg_e3[j];
            if (P[i][j] > 0 && indep > 0) {
                D_KL += P[i][j] * log2(P[i][j] / indep);
            }
        }
    }

    return D_KL;
}

// ==================== Main ====================
int main(void) {
    uint8_t alpha4 = 0xc;
    uint8_t beta4  = 0x8;
    uint8_t gamma8 = 0x41;
    uint8_t l8     = 0x8;
    uint8_t n8     = 0x61;
    uint8_t h8     = 0x41;

    double max_kl = -1.0;
    uint8_t max_xi = 0;

    printf("Computing KL divergence for all xi (0..255)...\n\n");

    for (int xi = 0; xi < 256; xi++) {
        uint8_t xi8 = (uint8_t)xi;
        double dkl = evaluate_independence_like_paper(alpha4, beta4, gamma8, xi8, l8, n8, h8);
        printf("xi = 0x%02X (%3d)  D_KL = %.8f\n", xi8, xi, dkl);
        if (dkl > max_kl) {
            max_kl = dkl;
            max_xi = xi8;
        }
    }

    printf("\n========== RESULT ==========\n");
    printf("Maximum KL divergence: D_KL = %.8f  (achieved at xi = 0x%02X)\n", max_kl, max_xi);
    return 0;
}