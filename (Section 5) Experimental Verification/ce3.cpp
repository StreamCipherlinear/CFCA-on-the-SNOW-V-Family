#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

/* Global lookup table: val[input1][input2][output] = automaton8 result (-100.0 if invalid).
   Size: 256*256*256 ~ 134 MB. */
double val[256][256][256];

double automaton8(int input1, int input2, int output) {
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

double automaton8_with_k(int input1, int input2, int output, int *k_out) {
    double k = 0;
    int e = 0, t, flag;
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
    *k_out = (int)k;
    return sign * pow(2.0, -k);
}

int main() {
    /* ===== Safe values for gamma, xi, h (excluding 0x00) ===== */
    uint8_t safe_vals[] = {
        0x17, 0x28, 0x3F, 0x41, 0x56, 0x69, 0x7E, 0x83, 0x94, 0xAB, 0xBC, 0xC2, 0xD5, 0xEA, 0xFD
    };
    int num_safe = sizeof(safe_vals) / sizeof(safe_vals[0]);  // 15

    int gamma, n, theta;
    int k;
    bool pair_found[256][256] = { false };

    /* Step 1: Find all (gamma, n) with k=1 */
    for (int i = 0; i < num_safe; i++) {
        gamma = safe_vals[i];
        for (n = 0x10; n <= 0xFF; n++) {
            for (theta = 0x10; theta <= 0xFF; theta++) {
                double res = automaton8_with_k(n, theta, gamma, &k);
                if (res == -100.0) continue;
                if (k == 1) {
                    pair_found[gamma][n] = true;
                    break;
                }
            }
        }
    }

    /* Step 2: Precompute lookup table */
    printf("Precomputing automaton8 lookup table...\n");
    for (int a = 0; a < 256; a++) {
        for (int b = 0; b < 256; b++) {
            for (int c = 0; c < 256; c++) {
                val[a][b][c] = automaton8(a, b, c);
            }
        }
    }
    printf("Lookup table ready.\n");

    /* Step 3: Global search for optimal (beta, gamma, xi, h, n) */
    double best_global_c = -1e300;
    int best_beta = 0, best_gamma = 0, best_xi = 0, best_h = 0, best_n = 0;

    double A[256];

    for (int i = 0; i < num_safe; i++) {
        gamma = safe_vals[i];
        for (n = 0x10; n <= 0xFF; n++) {
            if (!pair_found[gamma][n]) continue;

            // Precompute A[theta] = val[theta][n][gamma]
            for (theta = 0; theta < 256; theta++) {
                double v = val[theta][n][gamma];
                A[theta] = (v == -100.0) ? 0.0 : v;
            }

            // Search for optimal beta, xi, h for this (gamma, n)
            double best_local_c = -1e300;
            int best_beta_loc = 0, best_xi_loc = 0, best_h_loc = 0;

            for (int beta = 0x01; beta <= 0x0F; beta++) {
                for (int j = 0; j < num_safe; j++) {
                    int xi = safe_vals[j];
                    int x = beta ^ xi;
                    for (int k_idx = 0; k_idx < num_safe; k_idx++) {
                        int h = safe_vals[k_idx];

                        double c = 0.0;
                        for (theta = 0; theta < 256; theta++) {
                            double term1 = A[theta];
                            if (term1 == 0.0) continue;
                            double term2 = val[x][h][theta];
                            if (term2 == -100.0) continue;
                            c += term1 * term2;
                        }

                        if (c > best_local_c) {
                            best_local_c = c;
                            best_beta_loc = beta;
                            best_xi_loc = xi;
                            best_h_loc = h;
                        }
                    }
                }
            }

            // Update global best
            if (best_local_c > best_global_c) {
                best_global_c = best_local_c;
                best_beta = best_beta_loc;
                best_gamma = gamma;
                best_xi = best_xi_loc;
                best_h = best_h_loc;
                best_n = n;
            }
        }
    }

    /* Step 4: Output the single best combination */
    printf("\nOptimal combination (maximum c):\n");
    double c_val = best_global_c;
    int sign = (c_val > 0) ? 1 : (c_val < 0) ? -1 : 0;
    double abs_c = fabs(c_val);
    double log2_abs = (abs_c > 0) ? log2(abs_c) : -INFINITY;

    printf("c = %+.6f  (sign=%+d, log2(|c|)=%.3f)\n", c_val, sign, log2_abs);
    printf("beta = 0x%02X\n", best_beta);
    printf("gamma = 0x%02X\n", best_gamma);
    printf("xi = 0x%02X\n", best_xi);
    printf("h = 0x%02X\n", best_h);
    printf("n = 0x%02X\n", best_n);

    return 0;
}