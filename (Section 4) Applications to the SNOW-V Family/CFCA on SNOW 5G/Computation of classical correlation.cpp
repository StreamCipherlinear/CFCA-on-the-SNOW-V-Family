#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* ---------- GF(2^8) ---------- */

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

/* ---------- Lookup tables ---------- */

static signed char sbox_sign[256][256];
static signed char sign_dot[256][256];
static int16_t    P[256][256][256];
static double     W_tab[256][256];

void build_tables() {
    /* sbox_sign[v][z] = (-1)^{v . S(z)} */
    for (int v = 0; v < 256; v++)
        for (int z = 0; z < 256; z++)
            sbox_sign[v][z] = DOT_TABLE[v][S_BOX[z]] ? -1 : 1;

    /* sign_dot[l][u] = (-1)^{l . u} */
    for (int l = 0; l < 256; l++)
        for (int u = 0; u < 256; u++)
            sign_dot[l][u] = DOT_TABLE[l][u] ? -1 : 1;

    /* P[v][l][w] = sum_y (-1)^{v.S(y) XOR l.((w-y)&0xFF)} */
    for (int v = 0; v < 256; v++)
        for (int l = 0; l < 256; l++)
            for (int w = 0; w < 256; w++) {
                int sum = 0;
                for (int y = 0; y < 256; y++)
                    sum += sbox_sign[v][y] * sign_dot[l][(w - y) & 0xFF];
                P[v][l][w] = (int16_t)sum;
            }

    /* W_tab[h][a] = (1/256) sum_x (-1)^{h.S(x) XOR a.x} */
    for (int h = 0; h < 256; h++)
        for (int a = 0; a < 256; a++) {
            int sum = 0;
            for (int x = 0; x < 256; x++)
                sum += sbox_sign[h][x] * sign_dot[a][x];
            W_tab[h][a] = (double)sum / 256.0;
        }
}

/* Q[v][alpha][l] = sum_w P[v][l][w] * (-1)^{alpha.w} */
static double Q_from_P(int v, int alpha, int l) {
    const int16_t *Plw = P[v][l];
    const signed char *sa = sign_dot[alpha];
    int32_t sum = 0;
    for (int w = 0; w < 256; w++)
        sum += (int32_t)Plw[w] * sa[w];
    return (double)sum;
}

/* ---------- automaton16 ---------- */

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

/* ---------- 4x4 matrix build + diagonal extraction ---------- */

static void build_matrix(const uint32_t mask[4], uint8_t mat[4][4]) {
    for (int k = 0; k < 4; k++) {
        int c = 3 - k;
        mat[0][c] = (mask[k]      ) & 0xFF;
        mat[1][c] = (mask[k] >>  8) & 0xFF;
        mat[2][c] = (mask[k] >> 16) & 0xFF;
        mat[3][c] = (mask[k] >> 24) & 0xFF;
    }
}

/* d=3 diagonal: (0,3),(1,0),(2,1),(3,2) */
static void extract_d3(const uint8_t mat[4][4], uint8_t out[4]) {
    out[0] = mat[0][3]; out[1] = mat[1][0];
    out[2] = mat[2][1]; out[3] = mat[3][2];
}

/* d=1 diagonal: (0,1),(1,2),(2,3),(3,0) */
static void extract_d1(const uint8_t mat[4][4], uint8_t out[4]) {
    out[0] = mat[0][1]; out[1] = mat[1][2];
    out[2] = mat[2][3]; out[3] = mat[3][0];
}

/* ---------- Global state ---------- */

#define ITERATIONS 3

static uint8_t alpha_H[4], l_H[4];
static uint8_t alpha_L[4], l_L[4];
static uint8_t v1_H_glob[4], v1_L_glob[4];

static int valid_H[4], valid_L[4];
static double  e3_H_all[256], e3_L_all[256];
static uint8_t res2_H_all[256][4], res2_L_all[256][4];

/* corr(h, v2; alpha, l) = W_h(alpha) * Q_v2(alpha, l) / 65536 */
static double corr_value(int h, int v2, int alpha, int l) {
    return W_tab[h][alpha] * Q_from_P(v2, alpha, l) / 65536.0;
}

/* Compute S_H and S_L with current masks */
static void compute_S(double* out_SH, double* out_SL) {
    double SH = 0, SL = 0;
    for (int xi = 0; xi < 256; xi++) {
        double prod = 1.0; int ok = 1;
        for (int p = 0; p < 4; p++) {
            if (!valid_H[p]) continue;
            int v2 = res2_H_all[xi][p];
            if (v2 == 0) { ok = 0; break; }
            prod *= corr_value(v1_H_glob[p], v2, alpha_H[p], l_H[p]);
        }
        if (ok) SH += prod * e3_H_all[xi];
    }
    for (int xi = 0; xi < 256; xi++) {
        double prod = 1.0; int ok = 1;
        for (int p = 0; p < 4; p++) {
            if (!valid_L[p]) continue;
            int v2 = res2_L_all[xi][p];
            if (v2 == 0) { ok = 0; break; }
            prod *= corr_value(v1_L_glob[p], v2, alpha_L[p], l_L[p]);
        }
        if (ok) SL += prod * e3_L_all[xi];
    }
    *out_SH = SH; *out_SL = SL;
}

/* ---------- main ---------- */

int main() {
    clock_t t_start = clock();

    /* Phase 1: build tables (silent) */
    compute_dot_table();
    build_tables();

    /* Phase 2: v1 from h */
    const uint16_t h_H = 0xe1a0, h_L = 0x7f40;
    allbox_optimized((uint32_t)h_H, v1_H_glob);
    allbox_optimized((uint32_t)h_L, v1_L_glob);

    for (int p = 0; p < 4; p++) {
        valid_H[p] = (v1_H_glob[p] != 0);
        valid_L[p] = (v1_L_glob[p] != 0);
    }
    printf("[Phase 2] v1_H = %02x %02x %02x %02x, v1_L = %02x %02x %02x %02x\n",
           v1_H_glob[0],v1_H_glob[1],v1_H_glob[2],v1_H_glob[3],
           v1_L_glob[0],v1_L_glob[1],v1_L_glob[2],v1_L_glob[3]);

    /* Phase 3: per-xi e3 and v2 (silent) */
    for (int xi = 0; xi < 256; xi++) {
        uint16_t xi_H = 0x8100 | xi;
        e3_H_all[xi] = cor_e3_H(xi_H);
        allbox_optimized((uint32_t)xi_H, res2_H_all[xi]);
    }
    for (int xi = 0; xi < 256; xi++) {
        uint16_t xi_L = 0x7f00 | xi;
        e3_L_all[xi] = cor_e3_L(xi_L);
        allbox_optimized((uint32_t)xi_L, res2_L_all[xi]);
    }

    /* Phase 4: initial masks (silent) */
    const uint32_t ALPHA_INIT[4] = {0x00000002u, 0x0d004100u, 0x00020041u, 0x01000000u};
    const uint32_t L_INIT[4]     = {0x00000002u, 0x0d006100u, 0x00020061u, 0x01000000u};

    uint8_t mat_a[4][4], mat_l[4][4];
    build_matrix(ALPHA_INIT, mat_a);
    build_matrix(L_INIT, mat_l);
    extract_d3(mat_a, alpha_H); extract_d3(mat_l, l_H);
    extract_d1(mat_a, alpha_L); extract_d1(mat_l, l_L);

    double SH0, SL0;
    compute_S(&SH0, &SL0);

    /* Phase 5: coordinate descent (silent) */
    for (int iter = 0; iter < ITERATIONS; iter++) {
        /* ==== H side ==== */
        for (int p = 0; p < 4; p++) {
            if (!valid_H[p]) continue;

            /* other_H[xi] = e3_H(xi) * prod_{q != p} corr_q */
            static double other_H[256];
            for (int xi = 0; xi < 256; xi++) {
                double prod = e3_H_all[xi];
                int ok = 1;
                for (int q = 0; q < 4; q++) {
                    if (q == p || !valid_H[q]) continue;
                    int v2 = res2_H_all[xi][q];
                    if (v2 == 0) { ok = 0; break; }
                    prod *= corr_value(v1_H_glob[q], v2, alpha_H[q], l_H[q]);
                }
                other_H[xi] = ok ? prod : 0.0;
            }

            /* S2[l][w] = sum_xi other_H(xi) * P[v2(xi)[p]][l][w] */
            static double S2[256][256];
            for (int l = 0; l < 256; l++)
                for (int w = 0; w < 256; w++)
                    S2[l][w] = 0.0;

            for (int xi = 0; xi < 256; xi++) {
                if (other_H[xi] == 0.0) continue;
                int v2 = res2_H_all[xi][p];
                if (v2 == 0) continue;
                double ov = other_H[xi];
                const int16_t *Pv2 = &P[v2][0][0];
                for (int l = 0; l < 256; l++) {
                    const int16_t *Prow = Pv2 + l * 256;
                    double *Srow = S2[l];
                    for (int w = 0; w < 256; w++)
                        Srow[w] += ov * (double)Prow[w];
                }
            }

            /* SL_now: current S_L (fixed while optimizing H) */
            double SL_now = 0;
            for (int xi = 0; xi < 256; xi++) {
                double prod = 1.0; int ok = 1;
                for (int q = 0; q < 4; q++) {
                    if (!valid_L[q]) continue;
                    int v2 = res2_L_all[xi][q];
                    if (v2 == 0) { ok = 0; break; }
                    prod *= corr_value(v1_L_glob[q], v2, alpha_L[q], l_L[q]);
                }
                if (ok) SL_now += prod * e3_L_all[xi];
            }

            /* search (alpha, l) */
            double best_abs = -1;
            int best_alpha = alpha_H[p], best_l = l_H[p];
            for (int alpha = 0; alpha < 256; alpha++) {
                double Wh = W_tab[v1_H_glob[p]][alpha] / 65536.0;
                if (Wh == 0.0) continue;
                const signed char *sa = sign_dot[alpha];
                for (int l = 0; l < 256; l++) {
                    const double *S2l = S2[l];
                    double sumY = 0;
                    for (int w = 0; w < 256; w++)
                        sumY += S2l[w] * sa[w];
                    double aS = fabs(Wh * sumY * SL_now);
                    if (aS > best_abs) {
                        best_abs = aS;
                        best_alpha = alpha;
                        best_l = l;
                    }
                }
            }
            alpha_H[p] = (uint8_t)best_alpha;
            l_H[p] = (uint8_t)best_l;
        }

        compute_S(&SH0, &SL0);

        /* ==== L side ==== */
        for (int p = 0; p < 4; p++) {
            if (!valid_L[p]) continue;

            static double other_L[256];
            for (int xi = 0; xi < 256; xi++) {
                double prod = e3_L_all[xi];
                int ok = 1;
                for (int q = 0; q < 4; q++) {
                    if (q == p || !valid_L[q]) continue;
                    int v2 = res2_L_all[xi][q];
                    if (v2 == 0) { ok = 0; break; }
                    prod *= corr_value(v1_L_glob[q], v2, alpha_L[q], l_L[q]);
                }
                other_L[xi] = ok ? prod : 0.0;
            }

            static double S2[256][256];
            for (int l = 0; l < 256; l++)
                for (int w = 0; w < 256; w++)
                    S2[l][w] = 0.0;

            for (int xi = 0; xi < 256; xi++) {
                if (other_L[xi] == 0.0) continue;
                int v2 = res2_L_all[xi][p];
                if (v2 == 0) continue;
                double ov = other_L[xi];
                const int16_t *Pv2 = &P[v2][0][0];
                for (int l = 0; l < 256; l++) {
                    const int16_t *Prow = Pv2 + l * 256;
                    double *Srow = S2[l];
                    for (int w = 0; w < 256; w++)
                        Srow[w] += ov * (double)Prow[w];
                }
            }

            double SH_now = 0;
            for (int xi = 0; xi < 256; xi++) {
                double prod = 1.0; int ok = 1;
                for (int q = 0; q < 4; q++) {
                    if (!valid_H[q]) continue;
                    int v2 = res2_H_all[xi][q];
                    if (v2 == 0) { ok = 0; break; }
                    prod *= corr_value(v1_H_glob[q], v2, alpha_H[q], l_H[q]);
                }
                if (ok) SH_now += prod * e3_H_all[xi];
            }

            double best_abs = -1;
            int best_alpha = alpha_L[p], best_l = l_L[p];
            for (int alpha = 0; alpha < 256; alpha++) {
                double Wh = W_tab[v1_L_glob[p]][alpha] / 65536.0;
                if (Wh == 0.0) continue;
                const signed char *sa = sign_dot[alpha];
                for (int l = 0; l < 256; l++) {
                    const double *S2l = S2[l];
                    double sumY = 0;
                    for (int w = 0; w < 256; w++)
                        sumY += S2l[w] * sa[w];
                    double aS = fabs(SH_now * Wh * sumY);
                    if (aS > best_abs) {
                        best_abs = aS;
                        best_alpha = alpha;
                        best_l = l;
                    }
                }
            }
            alpha_L[p] = (uint8_t)best_alpha;
            l_L[p] = (uint8_t)best_l;
        }

        compute_S(&SH0, &SL0);
    }

    /* Phase 6: output log2|S| and final masks */
    double SH_final, SL_final;
    compute_S(&SH_final, &SL_final);
    double S_total = SH_final * SL_final;

    printf("[Phase 6] log2|S| : %.6f\n", log2(fabs(S_total)));
    printf("          alpha_H = %02x %02x %02x %02x, l_H = %02x %02x %02x %02x\n",
           alpha_H[0],alpha_H[1],alpha_H[2],alpha_H[3],
           l_H[0],l_H[1],l_H[2],l_H[3]);
    printf("          alpha_L = %02x %02x %02x %02x, l_L = %02x %02x %02x %02x\n",
           alpha_L[0],alpha_L[1],alpha_L[2],alpha_L[3],
           l_L[0],l_L[1],l_L[2],l_L[3]);

    /* Phase 7: optimal (xi^H, xi^L) */
    static double term_H[256], term_L[256];
    for (int xi = 0; xi < 256; xi++) {
        double prod = 1.0; int ok = 1;
        for (int p = 0; p < 4; p++) {
            if (!valid_H[p]) continue;
            int v2 = res2_H_all[xi][p];
            if (v2 == 0) { ok = 0; break; }
            prod *= corr_value(v1_H_glob[p], v2, alpha_H[p], l_H[p]);
        }
        term_H[xi] = ok ? prod * e3_H_all[xi] : 0.0;
    }
    for (int xi = 0; xi < 256; xi++) {
        double prod = 1.0; int ok = 1;
        for (int p = 0; p < 4; p++) {
            if (!valid_L[p]) continue;
            int v2 = res2_L_all[xi][p];
            if (v2 == 0) { ok = 0; break; }
            prod *= corr_value(v1_L_glob[p], v2, alpha_L[p], l_L[p]);
        }
        term_L[xi] = ok ? prod * e3_L_all[xi] : 0.0;
    }

    double best_abs_term = -1.0;
    int best_iH = -1, best_iL = -1;
    for (int i = 0; i < 256; i++) {
        if (term_H[i] == 0.0) continue;
        for (int j = 0; j < 256; j++) {
            if (term_L[j] == 0.0) continue;
            double aT = fabs(term_H[i] * term_L[j]);
            if (aT > best_abs_term) {
                best_abs_term = aT;
                best_iH = i;
                best_iL = j;
            }
        }
    }
    if (best_iH >= 0 && best_iL >= 0) {
        printf("[Phase 7] Optimal xi^H = 0x%04x, xi^L = 0x%04x\n",
               0x8100 | best_iH, 0x7f00 | best_iL);
        printf("          v2_H = %02x %02x %02x %02x\n",
               res2_H_all[best_iH][0], res2_H_all[best_iH][1],
               res2_H_all[best_iH][2], res2_H_all[best_iH][3]);
        printf("          v2_L = %02x %02x %02x %02x\n",
               res2_L_all[best_iL][0], res2_L_all[best_iL][1],
               res2_L_all[best_iL][2], res2_L_all[best_iL][3]);
        printf("          log2|term_xi| : %.6f\n", log2(best_abs_term));
    }

    /* Phase 8: only log2|term| at fixed xi = (0x8180, 0x7f40) */
    const uint32_t ALPHA_9[4] = {0x00000002u, 0x80006000u, 0x00020041u, 0xb0000000u};
    const uint32_t L_9[4]     = {0x00000002u, 0xc0006000u, 0x00020061u, 0xf0000000u};

    uint8_t mat_a9[4][4], mat_l9[4][4];
    build_matrix(ALPHA_9, mat_a9);
    build_matrix(L_9, mat_l9);

    uint8_t aH9[4], lH9[4], aL9[4], lL9[4];
    extract_d3(mat_a9, aH9); extract_d3(mat_l9, lH9);
    extract_d1(mat_a9, aL9); extract_d1(mat_l9, lL9);

    for (int p = 0; p < 4; p++) {
        alpha_H[p] = aH9[p]; l_H[p] = lH9[p];
        alpha_L[p] = aL9[p]; l_L[p] = lL9[p];
    }

    int iH0 = 0x8180 & 0xFF;
    int iL0 = 0x7f40 & 0xFF;

    double e2_H = 1.0;
    for (int p = 0; p < 4; p++) {
        if (!valid_H[p]) continue;
        int v2 = res2_H_all[iH0][p];
        if (v2 == 0) { e2_H = 0; break; }
        e2_H *= corr_value(v1_H_glob[p], v2, alpha_H[p], l_H[p]);
    }
    double e2_L = 1.0;
    for (int p = 0; p < 4; p++) {
        if (!valid_L[p]) continue;
        int v2 = res2_L_all[iL0][p];
        if (v2 == 0) { e2_L = 0; break; }
        e2_L *= corr_value(v1_L_glob[p], v2, alpha_L[p], l_L[p]);
    }
    double term = e2_H * e2_L * e3_H_all[iH0] * e3_L_all[iL0];

    printf("[Phase 8] log2|term| : %.6f\n", log2(fabs(term)));

    clock_t t_end = clock();
    printf("          Time: %.2f sec\n", (double)(t_end - t_start) / CLOCKS_PER_SEC);
    return 0;
}