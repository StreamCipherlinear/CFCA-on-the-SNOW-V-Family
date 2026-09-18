#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* ==================== GF(2^8) helpers ==================== */

uint8_t gfmul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (int i = 0; i < 8; i++) {
        if (b & 1) p ^= a;
        uint8_t hi = a & 0x80;
        a <<= 1;
        if (hi) a ^= 0x1b;
        b >>= 1;
    }
    return p;
}

uint8_t dot_product(uint8_t a, uint8_t b) {
    uint8_t r = a & b;
    r ^= (r >> 4); r ^= (r >> 2); r ^= (r >> 1);
    return r & 1;
}

static const uint8_t S_BOX[256] = {
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};

/* ==================== DOT table ==================== */

uint8_t DOT_TABLE[256][256];

void compute_dot_table() {
    for (int i = 0; i < 256; i++)
        for (int j = 0; j < 256; j++) {
            uint8_t r = 0;
            for (int k = 0; k < 8; k++)
                if (((i >> k) & 1) && ((j >> k) & 1)) r ^= 1;
            DOT_TABLE[i][j] = r;
        }
}

/* ==================== allbox_optimized ==================== */

void allbox_optimized(uint32_t t, uint8_t* results) {
    static uint8_t x_times_2[256];
    static int initialized = 0;
    if (!initialized) {
        for (int x = 0; x < 256; x++) x_times_2[x] = gfmul(x, 2);
        initialized = 1;
    }
    uint8_t t1 = t & 0xFF, t2 = (t >> 8) & 0xFF, t3 = (t >> 16) & 0xFF, t4 = (t >> 24) & 0xFF;
    uint8_t c1 = t1 ^ t4, c2 = t2 ^ t3 ^ t4, c3 = t1 ^ t2, c4 = t1 ^ t3 ^ t4,
            c5 = t2 ^ t3, c6 = t1 ^ t2 ^ t4, c7 = t3 ^ t4, c8 = t1 ^ t2 ^ t3;
    uint8_t cs[4][2] = {{c1,c2},{c3,c4},{c5,c6},{c7,c8}};
    for (int case_num = 0; case_num < 4; case_num++) {
        uint8_t ca = cs[case_num][0], cb = cs[case_num][1];
        uint8_t a = 0;
        for (int k = 0; k < 8; k++) {
            int x = 1 << k;
            int fx = dot_product(ca, x_times_2[x]) ^ dot_product(cb, x);
            if (fx) a |= (1 << k);
        }
        int linear = 1;
        for (int x = 0; x < 256; x++) {
            int fx = dot_product(ca, x_times_2[x]) ^ dot_product(cb, x);
            if (fx != dot_product(a, x)) { linear = 0; break; }
        }
        results[case_num] = linear ? a : 0;
    }
}

/* ==================== Precomputed tables ==================== */

static signed char sbox_sign[256][256];   /* (-1)^{v.S(z)} */
static signed char sign_dot[256][256];    /* (-1)^{l.u}    */
static int16_t    P[256][256][256];       /* P[v][l][w]    */

static int xi_list_H[4][256][256];
static int xi_list_H_cnt[4][256];
static int xi_list_L[4][256][256];
static int xi_list_L_cnt[4][256];

/* Per-s optimal l, saved from Step 6 */
static uint8_t lH_final[4][256];   /* [p][s] */
static uint8_t lL_final[4][256];

/* ==================== automaton16 ==================== */

double automaton16(uint32_t input1, uint32_t input2, uint32_t output) {
    double k = 0; int e = 0; double sign = 1.0;
    for (int t = 0; t < 16; t++) {
        int flag = ((output >> (15 - t)) & 1) * 4
                 + ((input1 >> (15 - t)) & 1) * 2
                 + ((input2 >> (15 - t)) & 1);
        if (e == 0) {
            if (flag == 7) e = 1;
            else if (flag == 0) e = 0;
            else return -100.0;
        } else {
            if (flag == 0 || flag == 3 || flag == 5 || flag == 6) { e = 1; k++; }
            else { e = 0; k++; }
            if (flag == 3 || flag == 4) sign = -sign;
        }
    }
    return sign * pow(2.0, -k);
}

double cor_e3_H(uint16_t xi_H) {
    double sum = 0.0;
    for (int i = 0; i < 256; i++) {
        uint32_t theta = i << 8;
        double p1 = automaton16(theta, 0x81e8, 0x81ec);
        double p2 = automaton16(xi_H, 0xe1a0, theta);
        if (p1 > -100 && p2 > -100) sum += p1 * p2;
    }
    return sum;
}

double cor_e3_L(uint16_t xi_L) {
    double sum = 0.0;
    for (int i = 0; i < 256; i++) {
        uint32_t theta = i << 8;
        double p1 = automaton16(theta, 0x5a00, 0x5a80);
        double p2 = automaton16(xi_L ^ 0x0040, 0x7f40, theta);
        if (p1 > -100 && p2 > -100) sum += p1 * p2;
    }
    return sum;
}

/* ==================== main ==================== */

#define COORD_ITER 3

int main() {
    clock_t t_start = clock();

    /* ---------- Step 1 ---------- */
    compute_dot_table();
    for (int l = 0; l < 256; l++)
        for (int u = 0; u < 256; u++)
            sign_dot[l][u] = DOT_TABLE[l][u] ? -1 : 1;
    for (int v = 0; v < 256; v++)
        for (int z = 0; z < 256; z++)
            sbox_sign[v][z] = DOT_TABLE[v][S_BOX[z]] ? -1 : 1;

    /* ---------- Step 2: P table ---------- */
    for (int v = 0; v < 256; v++)
        for (int l = 0; l < 256; l++)
            for (int w = 0; w < 256; w++) {
                int sum = 0;
                for (int x = 0; x < 256; x++)
                    sum += sbox_sign[v][x] * sign_dot[l][(w - x) & 0xFF];
                P[v][l][w] = (int16_t)sum;
            }

    /* Fixed parameters */
    const uint16_t h_H   = 0xe1a0;
    const uint16_t h_L   = 0x7f40;
    const uint16_t xi0_H = 0x8180;
    const uint16_t xi0_L = 0x7f40;

    /* ---------- Step 3: e3 and v2 tables ---------- */
    static double  e3_H_all[256], e3_L_all[256];
    static uint8_t res2_H_all[256][4], res2_L_all[256][4];

    for (int xi = 0; xi < 256; xi++) {
        uint16_t xi_H = 0x8100 | xi;
        e3_H_all[xi] = cor_e3_H(xi_H);
        uint8_t r2[4];
        allbox_optimized((uint32_t)xi_H, r2);
        for (int p = 0; p < 4; p++) res2_H_all[xi][p] = r2[p];
    }
    for (int xi = 0; xi < 256; xi++) {
        uint16_t xi_L = 0x7f00 | xi;
        e3_L_all[xi] = cor_e3_L(xi_L);
        uint8_t r2[4];
        allbox_optimized((uint32_t)xi_L, r2);
        for (int p = 0; p < 4; p++) res2_L_all[xi][p] = r2[p];
    }

    /* ---------- Step 4: fixed v1 ---------- */
    uint8_t v1_H[4], v1_L[4];
    allbox_optimized((uint32_t)h_H, v1_H);
    allbox_optimized((uint32_t)h_L, v1_L);

    uint8_t v2_H_0[4], v2_L_0[4];
    allbox_optimized((uint32_t)xi0_H, v2_H_0);
    allbox_optimized((uint32_t)xi0_L, v2_L_0);

    int valid_H[4], valid_L[4];
    for (int p = 0; p < 4; p++) {
        valid_H[p] = (v1_H[p] != 0 && v2_H_0[p] != 0);
        valid_L[p] = (v1_L[p] != 0 && v2_L_0[p] != 0);
    }

    /* ---------- Step 5: group xi by v2 ---------- */
    for (int p = 0; p < 4; p++) {
        for (int v2 = 0; v2 < 256; v2++) xi_list_H_cnt[p][v2] = 0;
        for (int xi = 0; xi < 256; xi++) {
            int v2 = res2_H_all[xi][p];
            xi_list_H[p][v2][xi_list_H_cnt[p][v2]++] = xi;
        }
    }
    for (int p = 0; p < 4; p++) {
        for (int v2 = 0; v2 < 256; v2++) xi_list_L_cnt[p][v2] = 0;
        for (int xi = 0; xi < 256; xi++) {
            int v2 = res2_L_all[xi][p];
            xi_list_L[p][v2][xi_list_L_cnt[p][v2]++] = xi;
        }
    }

    /* Scratch tables */
    static float corr_H[4][256][256];   /* [p][xi][l] */
    static float corr_L[4][256][256];

    /* Results */
    static double S_H_opt[256], S_L_opt[256];
    int tau_opt[256];

    /* ---------- Step 6: per-s optimization ---------- */
    int xi0_idx_H = xi0_H & 0xFF;
    int xi0_idx_L = xi0_L & 0xFF;

    for (int s = 0; s < 256; s++) {
        /* 6a: corr tables for this s */
        for (int p = 0; p < 4; p++) {
            if (!valid_H[p]) continue;
            int v1 = v1_H[p];
            const signed char *sv1 = sbox_sign[v1];
            int8_t sv1_s[256];
            for (int w = 0; w < 256; w++) sv1_s[w] = sv1[s ^ w];

            for (int v2 = 1; v2 < 256; v2++) {
                int cnt = xi_list_H_cnt[p][v2];
                if (cnt == 0) continue;
                int32_t accum[256];
                for (int l = 0; l < 256; l++) {
                    const int16_t *P_l = P[v2][l];
                    int32_t acc = 0;
                    for (int w = 0; w < 256; w++)
                        acc += (int32_t)P_l[w] * sv1_s[w];
                    accum[l] = acc;
                }
                for (int k = 0; k < cnt; k++) {
                    int xi = xi_list_H[p][v2][k];
                    float *out = corr_H[p][xi];
                    for (int l = 0; l < 256; l++)
                        out[l] = accum[l] / 65536.0f;
                }
            }
            for (int k = 0; k < xi_list_H_cnt[p][0]; k++) {
                int xi = xi_list_H[p][0][k];
                float *out = corr_H[p][xi];
                for (int l = 0; l < 256; l++) out[l] = 0.0f;
            }
        }
        for (int p = 0; p < 4; p++) {
            if (!valid_L[p]) continue;
            int v1 = v1_L[p];
            const signed char *sv1 = sbox_sign[v1];
            int8_t sv1_s[256];
            for (int w = 0; w < 256; w++) sv1_s[w] = sv1[s ^ w];

            for (int v2 = 1; v2 < 256; v2++) {
                int cnt = xi_list_L_cnt[p][v2];
                if (cnt == 0) continue;
                int32_t accum[256];
                for (int l = 0; l < 256; l++) {
                    const int16_t *P_l = P[v2][l];
                    int32_t acc = 0;
                    for (int w = 0; w < 256; w++)
                        acc += (int32_t)P_l[w] * sv1_s[w];
                    accum[l] = acc;
                }
                for (int k = 0; k < cnt; k++) {
                    int xi = xi_list_L[p][v2][k];
                    float *out = corr_L[p][xi];
                    for (int l = 0; l < 256; l++)
                        out[l] = accum[l] / 65536.0f;
                }
            }
            for (int k = 0; k < xi_list_L_cnt[p][0]; k++) {
                int xi = xi_list_L[p][0][k];
                float *out = corr_L[p][xi];
                for (int l = 0; l < 256; l++) out[l] = 0.0f;
            }
        }

        /* 6b: init l at xi_0 */
        uint8_t lH[4], lL[4];
        for (int p = 0; p < 4; p++) {
            if (!valid_H[p]) { lH[p] = 0; continue; }
            const float *c = corr_H[p][xi0_idx_H];
            float best = 0; int best_l = 0;
            for (int l = 0; l < 256; l++) {
                float ac = fabsf(c[l]);
                if (ac > best) { best = ac; best_l = l; }
            }
            lH[p] = (uint8_t)best_l;
        }
        for (int p = 0; p < 4; p++) {
            if (!valid_L[p]) { lL[p] = 0; continue; }
            const float *c = corr_L[p][xi0_idx_L];
            float best = 0; int best_l = 0;
            for (int l = 0; l < 256; l++) {
                float ac = fabsf(c[l]);
                if (ac > best) { best = ac; best_l = l; }
            }
            lL[p] = (uint8_t)best_l;
        }

        /* 6c: H descent */
        for (int iter = 0; iter < COORD_ITER; iter++) {
            for (int p = 0; p < 4; p++) {
                if (!valid_H[p]) continue;
                double best_abs = -1;
                int best_l = lH[p];
                for (int l = 0; l < 256; l++) {
                    double sum = 0;
                    for (int xi = 0; xi < 256; xi++) {
                        double prod = 1.0;
                        for (int q = 0; q < 4; q++) {
                            if (!valid_H[q]) continue;
                            int l_q = (q == p) ? l : lH[q];
                            prod *= corr_H[q][xi][l_q];
                        }
                        sum += prod * e3_H_all[xi];
                    }
                    double ac = fabs(sum);
                    if (ac > best_abs) { best_abs = ac; best_l = l; }
                }
                lH[p] = (uint8_t)best_l;
            }
        }

        /* 6d: L descent */
        for (int iter = 0; iter < COORD_ITER; iter++) {
            for (int p = 0; p < 4; p++) {
                if (!valid_L[p]) continue;
                double best_abs = -1;
                int best_l = lL[p];
                for (int l = 0; l < 256; l++) {
                    double sum = 0;
                    for (int xi = 0; xi < 256; xi++) {
                        double prod = 1.0;
                        for (int q = 0; q < 4; q++) {
                            if (!valid_L[q]) continue;
                            int l_q = (q == p) ? l : lL[q];
                            prod *= corr_L[q][xi][l_q];
                        }
                        sum += prod * e3_L_all[xi];
                    }
                    double ac = fabs(sum);
                    if (ac > best_abs) { best_abs = ac; best_l = l; }
                }
                lL[p] = (uint8_t)best_l;
            }
        }

        /* 6e: S_H[s], S_L[s] */
        double SH_s = 0, SL_s = 0;
        for (int xi = 0; xi < 256; xi++) {
            double prod = 1.0;
            for (int p = 0; p < 4; p++) {
                if (!valid_H[p]) continue;
                prod *= corr_H[p][xi][lH[p]];
            }
            SH_s += prod * e3_H_all[xi];
        }
        for (int xi = 0; xi < 256; xi++) {
            double prod = 1.0;
            for (int p = 0; p < 4; p++) {
                if (!valid_L[p]) continue;
                prod *= corr_L[p][xi][lL[p]];
            }
            SL_s += prod * e3_L_all[xi];
        }

        /* 6f: tau_s */
        tau_opt[s] = (SH_s * SL_s >= 0) ? 0 : 1;

        /* 6g: store */
        for (int p = 0; p < 4; p++) {
            lH_final[p][s] = lH[p];
            lL_final[p][s] = lL[p];
        }
        S_H_opt[s] = SH_s;
        S_L_opt[s] = SL_s;
    }

    /* ---------- Step 7: overall RMS ---------- */
    double sum_sq = 0;
    for (int s = 0; s < 256; s++) {
        double sign_tau = tau_opt[s] ? -1.0 : 1.0;
        double S = sign_tau * S_H_opt[s] * S_L_opt[s];
        sum_sq += S * S;
    }
    double rms = sqrt(sum_sq / 256.0);

    /* ---------- Step 8: find (xi^H, xi^L) maximizing RMS over s ---------- */
    static double term_H_all[256][256];   /* [s][xi] */
    static double term_L_all[256][256];

    for (int s = 0; s < 256; s++) {
        int8_t sv1_s_H[4][256];
        int8_t sv1_s_L[4][256];
        for (int p = 0; p < 4; p++) {
            if (valid_H[p])
                for (int w = 0; w < 256; w++)
                    sv1_s_H[p][w] = sbox_sign[v1_H[p]][s ^ w];
            if (valid_L[p])
                for (int w = 0; w < 256; w++)
                    sv1_s_L[p][w] = sbox_sign[v1_L[p]][s ^ w];
        }

        for (int xi = 0; xi < 256; xi++) {
            /* term_H[s][xi] */
            double prodH = 1.0; int okH = 1;
            for (int p = 0; p < 4; p++) {
                if (!valid_H[p]) continue;
                int v2 = res2_H_all[xi][p];
                if (v2 == 0) { okH = 0; break; }
                int l = lH_final[p][s];
                const int16_t *Pv2l = P[v2][l];
                int32_t sum = 0;
                for (int w = 0; w < 256; w++)
                    sum += (int32_t)Pv2l[w] * sv1_s_H[p][w];
                prodH *= sum / 65536.0;
            }
            term_H_all[s][xi] = okH ? prodH * e3_H_all[xi] : 0.0;

            /* term_L[s][xi] */
            double prodL = 1.0; int okL = 1;
            for (int p = 0; p < 4; p++) {
                if (!valid_L[p]) continue;
                int v2 = res2_L_all[xi][p];
                if (v2 == 0) { okL = 0; break; }
                int l = lL_final[p][s];
                const int16_t *Pv2l = P[v2][l];
                int32_t sum = 0;
                for (int w = 0; w < 256; w++)
                    sum += (int32_t)Pv2l[w] * sv1_s_L[p][w];
                prodL *= sum / 65536.0;
            }
            term_L_all[s][xi] = okL ? prodL * e3_L_all[xi] : 0.0;
        }
    }

    /* 8f: find (iH, iL) maximizing RMS over s */
    double best_rms = -1.0;
    int best_iH = -1, best_iL = -1;

    for (int iH = 0; iH < 256; iH++) {
        int any = 0;
        for (int s = 0; s < 256; s++)
            if (term_H_all[s][iH] != 0.0) { any = 1; break; }
        if (!any) continue;

        for (int iL = 0; iL < 256; iL++) {
            any = 0;
            for (int s = 0; s < 256; s++)
                if (term_L_all[s][iL] != 0.0) { any = 1; break; }
            if (!any) continue;

            double sum_sq_xi = 0;
            for (int s = 0; s < 256; s++) {
                double t = term_H_all[s][iH] * term_L_all[s][iL];
                sum_sq_xi += t * t;
            }
            double rms_xi = sqrt(sum_sq_xi / 256.0);
            if (rms_xi > best_rms) {
                best_rms = rms_xi;
                best_iH = iH;
                best_iL = iL;
            }
        }
    }

    clock_t t_end = clock();

    /* ---------- Output ---------- */
    printf("log2|RMS|       : %.6f\n", log2(fabs(rms)));

    if (best_iH >= 0 && best_iL >= 0) {
        printf("optimal xi^H    : 0x%04x\n", 0x8100 | best_iH);
        printf("optimal xi^L    : 0x%04x\n", 0x7f00 | best_iL);
        printf("v2_H(optimal)   : %02x %02x %02x %02x\n",
               res2_H_all[best_iH][0], res2_H_all[best_iH][1],
               res2_H_all[best_iH][2], res2_H_all[best_iH][3]);
        printf("v2_L(optimal)   : %02x %02x %02x %02x\n",
               res2_L_all[best_iL][0], res2_L_all[best_iL][1],
               res2_L_all[best_iL][2], res2_L_all[best_iL][3]);
        printf("log2|RMS_xi|    : %.6f\n", log2(best_rms));
    } else {
        printf("No valid (xi^H, xi^L) found.\n");
    }
    printf("Time            : %.2f sec\n", (double)(t_end - t_start) / CLOCKS_PER_SEC);

    return 0;
}