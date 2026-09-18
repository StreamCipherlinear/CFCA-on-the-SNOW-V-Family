#include <stdio.h>
#include <stdint.h>
#include <math.h>    
#include <stdlib.h>

// ========== Dot product ==========
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

int main() {
    // Fixed masks 
    uint8_t h  = 0x41;   
    uint8_t xi = 0x69;   

    printf("=== Adaptive l and tau for each u ===\n");
    printf("Fixed masks: h=0x%02X, xi=0x%02X\n", h, xi);
    printf("Expression: xi*R2 ^ l*y ^ h*R3  with R3 = u  ^ (x+y mod 16)\n");
    printf(" u  | best_l | tau |   cor\n");
    printf("----+--------+-----+----------\n");

    double sum_sq = 0.0;   // Sum of squared absolute correlations

    for (int u = 0; u < 16; u++) {
        int best_l = 0;
        double best_abs = 0.0;
        double best_cor = 0.0;

        for (int l = 0; l < 16; l++) {
            int s = 0;
            for (int x = 0; x < 16; x++) {
                for (int y = 0; y < 16; y++) {
                    int z = u ^ ((x + y) & 0xF);
                    if ((dot(xi, S8_transform(x)) ^ dot(l, y) ^ dot(h, S8_transform(z))) == 0)
                        s++;
                    else
                        s--;
                }
            }
            double cor = s / 256.0;
            if (fabs(cor) > best_abs) {
                best_abs = fabs(cor);
                best_cor = cor;
                best_l = l;
            }
        }

        int tau = (best_cor >= 0) ? 1 : 0;
        printf(" %2d |  0x%02X  |  %d  | %+9.6f\n", u, best_l, tau, best_cor);
        sum_sq += best_abs * best_abs;   // accumulate squared absolute correlation
    }

    // Compute root mean square (RMS) and its logarithm
    double rms = sqrt(sum_sq / 16.0);
    printf("\nRMS of |cor| over u = %.6f\n", rms);
    printf("log2(RMS) = %.6f\n", log2(rms));

    return 0;
}