/*
 * Search for optimal (n^L, xi^L, h^L) using cor(e1)*cor(e3^L).
 *
 * This program selects mask candidates for the lower 16-bit (+) block of e3
 * in SNOW 5G.  It combines the pre-computed correlation of e1 (read from
 * c(e1).txt) with the automaton output for P1^L and P2^L.
 *
 * The automaton evaluates:
 *   cor(P1^L)*cor(P2^L) = Sum_theta automaton(theta, n, gamma) * automaton(xi XOR beta, h, theta)
 *
 * For a fixed gamma = 0x5a80 (lower block), the program exhausts all combinations
 * of (beta, m) taken from c(e1).txt and all candidate masks from the predefined
 * sets:
 *   n  in {0x5200, 0x5300, 0x5a00, 0x5b00, 0x7200, 0x7300, 0x7a00, 0x7b00}
 *   xi, h  from value[table_idx][*]
 * The mask xi is XORed with beta before entering the automaton, reflecting the
 * structure (xi XOR beta) in e3^L.
 *
 * The intermediate mask theta is enumerated over its active byte (the lower byte,
 * i = 0x40..0x7f, with the upper byte zero).  Results with |cor| >= 2^-16 are
 * printed together with log2|cor| and sign and written to 0x5a80.txt.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

/*
 * Precomputed candidate sets for xi and h (lower block).
 * Each table contains 64 values corresponding to one active-byte pattern.
 * The four tables cover the possible patterns of the lower byte.
 */
uint32_t value[4][64] = {
    // value[0]
    {0x4080, 0x4183, 0x4284, 0x4387, 0x4489, 0x458a, 0x468d, 0x478e, 0x4891, 0x4992,
     0x4a95, 0x4b96, 0x4c98, 0x4d9b, 0x4e9c, 0x4f9f,
     0x50a0, 0x51a3, 0x52a4, 0x53a7, 0x54a9, 0x55aa, 0x56ad, 0x57ae, 0x58b1, 0x59b2,
     0x5ab5, 0x5bb6, 0x5cb8, 0x5dbb, 0x5ebc, 0x5fbf,
     0x60c0, 0x61c3, 0x62c4, 0x63c7, 0x64c9, 0x65ca, 0x66cd, 0x67ce, 0x68d1, 0x69d2,
     0x6ad5, 0x6bd6, 0x6cd8, 0x6ddb, 0x6edc, 0x6fdf,
     0x70e0, 0x71e3, 0x72e4, 0x73e7, 0x74e9, 0x75ea, 0x76ed, 0x77ee, 0x78f1, 0x79f2,
     0x7af5, 0x7bf6, 0x7cf8, 0x7dfb, 0x7efc, 0x7fff},
    // value[1]
    {0x403f, 0x41c0, 0x423e, 0x43c1, 0x443c, 0x45c3, 0x463d, 0x47c2, 0x48c7, 0x4938,
     0x4ac6, 0x4b39, 0x4cc4, 0x4d3b, 0x4ec5, 0x4f3a,
     0x5030, 0x51cf, 0x5231, 0x53ce, 0x5433, 0x55cc, 0x5632, 0x57cd, 0x58c8, 0x5937,
     0x5ac9, 0x5b36, 0x5ccb, 0x5d34, 0x5eca, 0x5f35,
     0x6020, 0x61df, 0x6221, 0x63de, 0x6423, 0x65dc, 0x6622, 0x67dd, 0x68d8, 0x6927,
     0x6ad9, 0x6b26, 0x6cdb, 0x6d24, 0x6eda, 0x6f25,
     0x702f, 0x71d0, 0x722e, 0x73d1, 0x742c, 0x75d3, 0x762d, 0x77d2, 0x78d7, 0x7928,
     0x7ad6, 0x7b29, 0x7cd4, 0x7d2b, 0x7ed5, 0x7f2a},
    // value[2]
    {0x4060, 0x41e1, 0x42e3, 0x4362, 0x4466, 0x45e7, 0x46e5, 0x4764, 0x48ec, 0x496d,
     0x4a6f, 0x4bee, 0x4cea, 0x4d6b, 0x4e69, 0x4fe8,
     0x50f8, 0x5179, 0x527b, 0x53fa, 0x54fe, 0x557f, 0x567d, 0x57fc, 0x5874, 0x59f5,
     0x5af7, 0x5b76, 0x5c72, 0x5df3, 0x5ef1, 0x5f70,
     0x6050, 0x61d1, 0x62d3, 0x6352, 0x6456, 0x65d7, 0x66d5, 0x6754, 0x68dc, 0x695d,
     0x6a5f, 0x6bde, 0x6cda, 0x6d5b, 0x6e59, 0x6fd8,
     0x70c8, 0x7149, 0x724b, 0x73ca, 0x74ce, 0x754f, 0x764d, 0x77cc, 0x7844, 0x79c5,
     0x7ac7, 0x7b46, 0x7c42, 0x7dc3, 0x7ec1, 0x7f40},
    // value[3]
    {0x4040, 0x4141, 0x4242, 0x4343, 0x4444, 0x4545, 0x4646, 0x4747, 0x4848, 0x4949,
     0x4a4a, 0x4b4b, 0x4c4c, 0x4d4d, 0x4e4e, 0x4f4f,
     0x5050, 0x5151, 0x5252, 0x5353, 0x5454, 0x5555, 0x5656, 0x5757, 0x5858, 0x5959,
     0x5a5a, 0x5b5b, 0x5c5c, 0x5d5d, 0x5e5e, 0x5f5f,
     0x6060, 0x6161, 0x6262, 0x6363, 0x6464, 0x6565, 0x6666, 0x6767, 0x6868, 0x6969,
     0x6a6a, 0x6b6b, 0x6c6c, 0x6d6d, 0x6e6e, 0x6f6f,
     0x7070, 0x7171, 0x7272, 0x7373, 0x7474, 0x7575, 0x7676, 0x7777, 0x7878, 0x7979,
     0x7a7a, 0x7b7b, 0x7c7c, 0x7d7d, 0x7e7e, 0x7f7f}
};

/**
 * 16-bit modular addition automaton.
 *
 * @param input1   first input mask  (v)
 * @param input2   second input mask (w)
 * @param output   output mask       (u)
 * @return signed correlation, or -100.0 if zero correlation is forced
 */
double automaton16(uint32_t input1, uint32_t input2, uint32_t output)
{
    double k = 0;
    int e = 0;
    int t;
    int flag;
    double sign = 1.0;

    for (t = 0; t < 16; t++) {
        flag = ((output >> (15 - t)) & 1) * 4 +
               ((input1 >> (15 - t)) & 1) * 2 +
               ((input2 >> (15 - t)) & 1);

        if (e == 0) {
            if (flag == 7)          e = 1;
            else if (flag == 0)     e = 0;
            else                    return -100.0;
        } else {
            if ((flag == 0) || (flag == 3) || (flag == 5) || (flag == 6))
                e = 1;
            else
                e = 0;
            k++;
            if ((flag == 3) || (flag == 4)) sign = -sign;
        }
    }
    return sign * pow(2.0, -k);
}

/* Candidate masks for n in the lower block (L) */
uint32_t N[8] = {0x5200, 0x5300, 0x5a00, 0x5b00,
                 0x7200, 0x7300, 0x7a00, 0x7b00};

int main()
{
    uint32_t gamma = 0x5a80;   // fixed output mask for the lower block

    /* Read pre-computed e1 correlations: beta, m, log2|cor|, sign */
    FILE *fp = fopen("c(e1).txt", "r");
    if (!fp) return 1;

    char filename[20];
    sprintf(filename, "0x%04x.txt", gamma);
    FILE *fp1 = fopen(filename, "w");

    int beta, m;
    double exponent;
    int sign;

    while (fscanf(fp, "%x %x %lf %d", &beta, &m, &exponent, &sign) == 4) {
        /* For each (beta, m) pair, search over mask candidates */
        #pragma omp parallel for collapse(4)
        for (int table_idx = 0; table_idx < 4; table_idx++)
        {
            for (int i4 = 0; i4 < 8; i4++)
            {
                for (int i3 = 0; i3 < 64; i3++)
                {
                    for (int i2 = 0; i2 < 64; i2++)
                    {
                        uint32_t xi = value[table_idx][i2];
                        uint32_t h  = value[table_idx][i3];
                        uint32_t n  = N[i4];
                        double ce3 = 0.0;

                        /* Enumerate the active low byte of theta (range 0x40..0x7F) */
                        #pragma omp simd reduction(+:ce3)
                        for (int i1 = 0x40; i1 < 128; i1++) {
                            uint32_t theta = i1 << 8;   // low byte active, high byte zero
                            double c1 = automaton16(theta, n, gamma);
                            double c2 = automaton16(xi ^ beta, h, theta);
                            if ((c1 > -100) && (c2 > -100))
                                ce3 += c1 * c2;
                        }

                        double cor = sign * pow(2.0, exponent) * ce3;
                        if (fabs(cor) >= pow(2.0, -16)) {
                            double cor_exp = log2(fabs(cor));
                            int cor_sign = (cor > 0) ? 1 : -1;
                            #pragma omp critical
                            {
                                printf("%d, 0x%04x, 0x%04x, 0x%04x, 0x%02x, 0x%02x, %.4f, %d\n",
                                       table_idx, n, xi, h, beta, m, cor_exp, cor_sign);
                                fprintf(fp1, "%d, 0x%04x, 0x%04x, 0x%04x, 0x%02x, 0x%02x, %.4f, %d\n",
                                        table_idx, n, xi, h, beta, m, cor_exp, cor_sign);
                            }
                        }
                    }
                }
            }
        }
    }

    fclose(fp);
    fclose(fp1);
    return 0;
}