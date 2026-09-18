#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

/* ========== From ce1e2.cpp ========== */

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

static double compute_e1_correlation(uint8_t beta, uint8_t m, uint8_t gamma) {
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

static double compute_e2_correlation(uint8_t alpha, uint8_t h, uint8_t xi, uint8_t l) {
    int s1 = 0;
    for (int x = 0; x < 16; x++) {
        if ((dot(alpha, x) ^ dot(h, S8_transform(x))) == 0)
            s1++;
        else
            s1--;
    }
    double cor1 = s1 / 16.0;

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
    return cor1 * cor2;
}

/* ========== From ce3.cpp ========== */

static double automaton8(int input1, int input2, int output) {
    double k = 0;
    int e = 0;
    int t, flag;
    double sign = 1.0;

    for (t = 0; t < 8; t++) {
        flag = ((output >> (7 - t)) & 0x1) * 4 +
               ((input1 >> (7 - t)) & 0x1) * 2 +
               ((input2 >> (7 - t)) & 0x1);

        if (e == 0) {
            if (flag == 7)      e = 1;
            else if (flag == 0) e = 0;
            else                return -100.0;
        } else {
            if (flag == 0 || flag == 3 || flag == 5 || flag == 6) {
                e = 1; k++;
            } else {
                e = 0; k++;
            }
            if (flag == 3 || flag == 4) sign = -sign;
        }
    }
    return sign * pow(2.0, -k);
}

static double compute_e3_correlation(uint8_t beta, uint8_t gamma, uint8_t h, uint8_t n, uint8_t xi) {
    double sum = 0.0;
    for (int theta = 0; theta < 256; theta++) {
        double aval = automaton8(theta, n, gamma);
        if (aval == -100.0) continue;
        double bval = automaton8(beta ^ xi, h, theta);
        if (bval == -100.0) continue;
        sum += aval * bval;
    }
    return sum;
}

/* ========== Main ========== */
int main() {
    const uint8_t BETA  = 0x08;
    const uint8_t M     = 0x08;
    const uint8_t GAMMA = 0x41;
    const uint8_t H     = 0x41;
    const uint8_t N     = 0x61;
    const uint8_t ALPHA = 0x0C;
    const uint8_t L     = 0x08;

    double cor_e1 = compute_e1_correlation(BETA, M, GAMMA);
    double sum_over_xi = 0.0;

    for (int xi = 0; xi < 256; xi++) {
        double cor_e2 = compute_e2_correlation(ALPHA, H, (uint8_t)xi, L);
        double cor_e3 = compute_e3_correlation(BETA, GAMMA, H, N, (uint8_t)xi);
        sum_over_xi += cor_e2 * cor_e3;
    }

    double cN = cor_e1 * sum_over_xi;
    double abs_cN = fabs(cN);
    double log2_abs = (abs_cN > 0.0) ? log2(abs_cN) : -INFINITY;

    printf("c(N) = %+.9f\n", cN);
    printf("log2(|c(N)|) = %.6f\n", log2_abs);

    return 0;
}