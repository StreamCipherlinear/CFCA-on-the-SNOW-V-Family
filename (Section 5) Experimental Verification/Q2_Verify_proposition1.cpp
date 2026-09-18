#include <stdio.h>
#include <stdint.h>
#include <math.h>

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

static void wht16_int(int *data) {
    for (int len = 1; len < 16; len <<= 1)
        for (int i = 0; i < 16; i += (len << 1))
            for (int j = 0; j < len; j++) {
                int a = data[i + j], b = data[i + j + len];
                data[i + j] = a + b;
                data[i + j + len] = a - b;
            }
}

int main() {
    uint8_t h = 0x41;
    int all_xi_pass = 1;

    for (int xi_int = 0; xi_int < 256; xi_int++) {
        uint8_t xi = (uint8_t)xi_int;
        int all_l_pass = 1;

        for (int l = 0; l < 16 && all_l_pass; l++) {
            int cor_int[16];

            for (int alpha = 0; alpha < 16; alpha++) {
                int s_alpha = 0;
                for (int x = 0; x < 16; x++)
                    for (int y = 0; y < 16; y++)
                        for (int u = 0; u < 16; u++) {
                            int z = u ^ ((x + y) & 0xF);
                            int g = dot(xi, S8_transform(x)) ^ dot(l, y) ^ dot(h, S8_transform(z));
                            int f = dot(alpha, u) ^ g;
                            s_alpha += (f == 0) ? 1 : -1;
                        }
                cor_int[alpha] = s_alpha;
            }

            wht16_int(cor_int);

            for (int u0 = 0; u0 < 16; u0++) {
                int s_left = 0;
                for (int x = 0; x < 16; x++)
                    for (int y = 0; y < 16; y++) {
                        int z = u0 ^ ((x + y) & 0xF);
                        int g = dot(xi, S8_transform(x)) ^ dot(l, y) ^ dot(h, S8_transform(z));
                        s_left += (g == 0) ? 1 : -1;
                    }
                if (s_left * 16 != cor_int[u0]) {
                    all_l_pass = 0;
                    break;
                }
            }
        }

        if (!all_l_pass) {
            all_xi_pass = 0;
            break;
        }
    }

    printf("%s\n", all_xi_pass ? "True" : "False");
    return 0;
}