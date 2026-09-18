/*
 * Correlation computation for e3 (SNOW-V family).
 *
 * This program evaluates the third sub-approximation e3 for given
 * linear mask candidates. It reads pre-computed e1 correlations
 * from "c(e1).txt" and exhaustively tests combinations of the masks
 * n, h, xi (together with beta, m) to filter candidates whose overall
 * correlation exceeds a minimum threshold.
 *
 * The automaton for 32-bit modular addition (+ mod 2^32) is used to compute
 * the correlation contributions. Look-up tables for row-vector
 * multiplications in the AES field GF(2^8) are precomputed to speed up
 * the generation of linear masks.
 *
 * Output is written to "0x81ec5a80.txt" (the filename reflects the
 * fixed output mask gamma).
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

/*
 * Precomputed candidate sets for the mask h (4 tables x 128 values).
 * Each table corresponds to a different active-byte pattern.
 * Values are 32-bit masks where the low 16 bits are zero.
 */
uint32_t value1[4][128] = {
    { /* table 0 */
     0x80c00000, 0x81410000, 0x82430000, 0x83c20000, 0x84c60000, 0x85470000, 0x86450000, 0x87c40000,
     0x884c0000, 0x89cd0000, 0x8acf0000, 0x8b4e0000, 0x8c4a0000, 0x8dcb0000, 0x8ec90000, 0x8f480000,
     /* ... remaining entries omitted for brevity */
     0xf0080000, 0xf1890000, 0xf28b0000, 0xf30a0000, 0xf40e0000, 0xf58f0000, 0xf68d0000, 0xf70c0000,
     0xf8840000, 0xf9050000, 0xfa070000, 0xfb860000, 0xfc820000, 0xfd030000, 0xfe010000, 0xff800000},
    { /* table 1 */
     0x80800000, 0x81810000, 0x82820000, 0x83830000, 0x84840000, 0x85850000, 0x86860000, 0x87870000,
     /* ... */
     0xf8f80000, 0xf9f90000, 0xfafa0000, 0xfbfb0000, 0xfcfc0000, 0xfdfd0000, 0xfefe0000, 0xffff0000},
    { /* table 2 */
     0x80010000, 0x81020000, 0x82050000, 0x83060000, 0x84080000, 0x850b0000, 0x860c0000, 0x870f0000,
     /* ... */
     0xf8f00000, 0xf9f30000, 0xfaf40000, 0xfbf70000, 0xfcf90000, 0xfdfa0000, 0xfefd0000, 0xfffe0000},
    { /* table 3 */
     0x807f0000, 0x81800000, 0x827e0000, 0x83810000, 0x847c0000, 0x85830000, 0x867d0000, 0x87820000,
     /* ... */
     0xf8a80000, 0xf9570000, 0xfaa90000, 0xfb560000, 0xfcab0000, 0xfd540000, 0xfeaa0000, 0xff550000}
};

/* Candidate set for the linear mask n (8 values). */
static const uint32_t N[8] = {
    0x91ec5a00, 0xa1ec5a00, 0xc1ec5a00, 0x81ec5a00,
    0xb1ec5a00, 0xd1ec5a00, 0xe1ec5a00, 0xf1ec5a00
};

/* Look-up tables for row-vector right-multiplication in GF(2^8). */
static uint8_t mul02_rowvector_table[256];
static uint8_t mul03_rowvector_table[256];
static uint8_t mul02inv_rowvector_table[256];
static uint8_t mul03inv_rowvector_table[256];

/*
 * Precompute row-vector right-multiplication tables for GF(2^8).
 * These tables encode the linear transformations that correspond to
 * multiplication by 2, 3, and their inverses when the input is
 * interpreted as a row vector (as opposed to the usual column vector
 * representation used in the AES specification).
 */
static void init_rowvector_tables(void)
{
    for (int i = 0; i < 256; i++) {
        uint8_t v = (uint8_t)i;

        uint8_t v7 = (v >> 7) & 1;
        uint8_t v6 = (v >> 6) & 1;
        uint8_t v5 = (v >> 5) & 1;
        uint8_t v4 = (v >> 4) & 1;
        uint8_t v3 = (v >> 3) & 1;
        uint8_t v2 = (v >> 2) & 1;
        uint8_t v1 = (v >> 1) & 1;
        uint8_t v0 = v & 1;

        /* ---- mul02 row-vector ---- */
        uint8_t c7 = (uint8_t)(v0 ^ v1 ^ v3 ^ v4);
        uint8_t c6 = v7;
        uint8_t c5 = v6;
        uint8_t c4 = v5;
        uint8_t c3 = v4;
        uint8_t c2 = v3;
        uint8_t c1 = v2;
        uint8_t c0 = v1;

        mul02_rowvector_table[i] = (uint8_t)((c7 << 7) | (c6 << 6) | (c5 << 5) | (c4 << 4) |
                                             (c3 << 3) | (c2 << 2) | (c1 << 1) | c0);
        /* mul03 = mul02 XOR identity */
        mul03_rowvector_table[i] = (uint8_t)(mul02_rowvector_table[i] ^ v);

        /* ---- mul02 inverse row-vector ---- */
        c7 = v6; c6 = v5; c5 = v4; c4 = v3; c3 = v2; c2 = v1; c1 = v0;
        c0 = (uint8_t)(v7 ^ v3 ^ v2 ^ v0);

        mul02inv_rowvector_table[i] = (uint8_t)((c7 << 7) | (c6 << 6) | (c5 << 5) | (c4 << 4) |
                                                (c3 << 3) | (c2 << 2) | (c1 << 1) | c0);

        /* ---- mul03 inverse row-vector ---- */
        c7 = (uint8_t)(v7 ^ v3 ^ v0);
        c6 = (uint8_t)(v7 ^ v6 ^ v3 ^ v0);
        c5 = (uint8_t)(v7 ^ v6 ^ v5 ^ v3 ^ v0);
        c4 = (uint8_t)(v7 ^ v6 ^ v5 ^ v4 ^ v3 ^ v0);
        c3 = (uint8_t)(v7 ^ v6 ^ v5 ^ v4 ^ v0);
        c2 = (uint8_t)(v7 ^ v6 ^ v5 ^ v4 ^ v2 ^ v0);
        c1 = (uint8_t)(v7 ^ v6 ^ v5 ^ v4 ^ v2 ^ v1 ^ v0);
        c0 = (uint8_t)(v7 ^ v6 ^ v5 ^ v4 ^ v2 ^ v1);

        mul03inv_rowvector_table[i] = (uint8_t)((c7 << 7) | (c6 << 6) | (c5 << 5) | (c4 << 4) |
                                                (c3 << 3) | (c2 << 2) | (c1 << 1) | c0);
    }
}

/* Inline look-up functions for the precomputed tables. */
static inline uint8_t mul02(uint8_t v)    { return mul02_rowvector_table[v]; }
static inline uint8_t mul03(uint8_t v)    { return mul03_rowvector_table[v]; }
static inline uint8_t mul02inv(uint8_t v) { return mul02inv_rowvector_table[v]; }
static inline uint8_t mul03inv(uint8_t v) { return mul03inv_rowvector_table[v]; }

/*
 * Compute the least significant byte of xi depending on the table index
 * and the three most significant bytes v3, v2, v1.
 *
 * The four cases correspond to different linear combinations of the
 * input bytes that arise from the structure of the mask template.
 */
static inline uint8_t xi_last_byte(int table_idx, uint8_t v3, uint8_t v2, uint8_t v1)
{
    switch (table_idx) {
        case 0: return mul02inv((uint8_t)(v1 ^ v2 ^ mul03(v3)));
        case 1: return mul03inv((uint8_t)(mul02(v1) ^ v2 ^ v3));
        case 2: return (uint8_t)(mul03(v1) ^ mul02(v2) ^ v3);
        case 3: return (uint8_t)(v1 ^ mul03(v2) ^ mul02(v3));
        default: return 0;
    }
}

/*
 * Assemble a 32-bit mask xi from three given bytes (v3, v2, v1) and
 * the table-dependent least significant byte.
 */
static inline uint32_t make_xi32(int table_idx, uint8_t v3, uint8_t v2, uint8_t v1)
{
    uint8_t x0 = xi_last_byte(table_idx, v3, v2, v1);
    return ((uint32_t)v3 << 24) | ((uint32_t)v2 << 16) | ((uint32_t)v1 << 8) | (uint32_t)x0;
}

/*
 * 32-bit modular addition automaton (+ mod 2^32).
 *
 * @param input1   first input mask  (v)
 * @param input2   second input mask (w)
 * @param output   output mask       (u)
 * @return signed correlation value, or -100.0 if zero correlation is forced.
 */
static inline double automaton32(uint32_t input1, uint32_t input2, uint32_t output)
{
    int k = 0;              /* exponent counter (solid-arrow transitions) */
    int e = 0;              /* automaton state (0 or 1) */
    double sign = 1.0;      /* sign factor */

    for (int t = 0; t < 32; t++) {
        int flag = ((output >> (31 - t)) & 1) * 4
                 + ((input1 >> (31 - t)) & 1) * 2
                 + ((input2 >> (31 - t)) & 1);

        if (e == 0) {
            if (flag == 7) e = 1;
            else if (flag == 0) e = 0;
            else return -100.0;          /* zero correlation forced */
        } else {
            k++;
            e = ((flag == 0) || (flag == 3) || (flag == 5) || (flag == 6)) ? 1 : 0;
            if ((flag == 3) || (flag == 4)) sign = -sign;
        }
    }
    return sign * ldexp(1.0, -k);         /* sign x 2^{-k} */
}

int main(void)
{
    init_rowvector_tables();

    /* Fixed output mask for the current search */
    const uint32_t gamma = 0x81ec5a80;
    const double THRESH = ldexp(1.0, -32);   /* minimum |cor| to report */

    /* Pre-build the list of theta values (active high byte, low three bytes zero) */
    uint32_t theta_list[128];
    for (int ti = 0; ti < 128; ti++) {
        theta_list[ti] = ((uint32_t)(ti + 0x80)) << 24;
    }

    /* Read pre-computed e1 correlations */
    FILE *fp = fopen("c(e1).txt", "r");
    if (!fp) { perror("fopen c(e1).txt"); return 1; }

    char filename[32];
    snprintf(filename, sizeof(filename), "0x%08x.txt", gamma);
    FILE *fp1 = fopen(filename, "w");
    if (!fp1) { perror("fopen output"); fclose(fp); return 1; }

    /* Pre-compute c1 = automaton(theta, n, gamma) for all n and theta */
    double c1_table[8][128];
    for (int i4 = 0; i4 < 8; i4++) {
        for (int ti = 0; ti < 128; ti++) {
            c1_table[i4][ti] = automaton32(theta_list[ti], N[i4], gamma);
        }
    }

    unsigned int beta, m;
    double exponent;
    int sign_in;

    /* Process each (beta, m) pair from the e1 correlation file */
    while (fscanf(fp, "%x %x %lf %d", &beta, &m, &exponent, &sign_in) == 4)
    {
        const uint32_t beta_u = (uint32_t)beta;
        const unsigned int m_u = m;
        const double scale = (double)sign_in * pow(2.0, exponent);

        /*
         * Parallel enumeration over table index, n, h, and the first byte of xi.
         * The inner two loops (over vv2, vv1) are kept serial to avoid
         * excessive parallel overhead.
         */
        #pragma omp parallel for collapse(4) schedule(static)
        for (int table_idx = 0; table_idx < 4; table_idx++)
            for (int i4 = 0; i4 < 8; i4++)
                for (int i3 = 0; i3 < 128; i3++)
                    for (int vv3 = 0x80; vv3 < 256; vv3++)
                    {
                        const uint32_t h = value1[table_idx][i3];

                        for (int vv2 = 0; vv2 < 256; vv2++)
                            for (int vv1 = 0; vv1 < 256; vv1++)
                            {
                                const uint32_t xi  = make_xi32(table_idx,
                                                               (uint8_t)vv3,
                                                               (uint8_t)vv2,
                                                               (uint8_t)vv1);
                                const uint32_t xib = xi ^ beta_u;

                                /* Accumulate c1 * c2 over all theta values */
                                double ce3_xi = 0.0;
                                for (int ti = 0; ti < 128; ti++) {
                                    const double c1 = c1_table[i4][ti];
                                    if (c1 <= -100.0) continue;

                                    const uint32_t theta = theta_list[ti];
                                    const double c2 = automaton32(xib, h, theta);
                                    if (c2 <= -100.0) continue;

                                    ce3_xi += c1 * c2;
                                }
                                if (ce3_xi == 0.0) continue;

                                const double cor = scale * ce3_xi;
                                if (fabs(cor) < THRESH) continue;

                                const double cor_exp  = log2(fabs(cor));
                                const int    cor_sign = (cor > 0.0) ? 1 : -1;

                                /* Output under mutual exclusion */
                                #pragma omp critical
                                {
                                    printf("%d, 0x%08x, 0x%08x, 0x%08x, 0x%02x, 0x%02x, %.4f, %d\n",
                                           table_idx, N[i4], h, xi, beta_u, m_u, cor_exp, cor_sign);
                                    fprintf(fp1, "%d, 0x%08x, 0x%08x, 0x%08x, 0x%02x, 0x%02x, %.4f, %d\n",
                                            table_idx, N[i4], h, xi, beta_u, m_u, cor_exp, cor_sign);
                                }
                            }
                    }
    }

    fclose(fp);
    fclose(fp1);
    return 0;
}