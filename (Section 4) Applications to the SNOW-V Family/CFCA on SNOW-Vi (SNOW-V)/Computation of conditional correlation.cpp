/*
 * Computation of the conditional correlation for SNOW-Vi (and SNOW-V).
 *
 * This program reads pre-screened mask candidates from the file
 * "0x81ec5a80.txt" (produced by the e3 filtering step) and evaluates the
 * conditional correlation of e2.  For each candidate tuple (h, xi) and each
 * active-byte position, it computes the eight best l-masks for every possible
 * keystream byte s (0..255).  The quality of a byte pair is summarised by the
 * log2 of the root mean square (RMS) of the maximal absolute correlations over
 * all s and over the top-ranked l values.  The program then combines the
 * per-byte exponents (8 choices per byte) and selects the combination that
 * yields the highest overall exponent (base exponent from the e1*e3 stage
 * plus the e2 contributions).  Only the globally best combination is reported.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

/* ------------------------- GF(2^8) arithmetic ----------------------- */

/**
 * Multiplication in GF(2^8) with the AES irreducible polynomial
 * x^8 + x^4 + x^3 + x + 1.
 */
uint8_t gfmul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    uint8_t hi_bit_set;

    for (int i = 0; i < 8; i++) {
        if (b & 1) p ^= a;
        hi_bit_set = a & 0x80;
        a <<= 1;
        if (hi_bit_set) a ^= 0x1b;   // reduction
        b >>= 1;
    }
    return p;
}

/**
 * Dot product (parity of bitwise AND).
 * Returns 0 if the number of common 1-bits is even, 1 otherwise.
 */
uint8_t dot_product(uint8_t a, uint8_t b) {
    uint8_t result = a & b;
    result ^= (result >> 4);
    result ^= (result >> 2);
    result ^= (result >> 1);
    return result & 0x01;
}

/**
 * For a given 32-bit mask t = (t4, t3, t2, t1) (big-endian byte order),
 * find the 8-bit output masks that yield a perfect correlation (+/-256)
 * for the four AES S-box equations (cases 0..3), if they exist.
 * Results are stored in results[0..3]; a zero value indicates that no
 * mask yields perfect correlation.
 */
void allbox_optimized(uint32_t t, uint8_t* results) {
    uint8_t t1 = t & 0xFF;
    uint8_t t2 = (t >> 8) & 0xFF;
    uint8_t t3 = (t >> 16) & 0xFF;
    uint8_t t4 = (t >> 24) & 0xFF;

    // pre-compute x*2 for all x
    uint8_t x_times_2[256];
    for (int x = 0; x < 256; x++) {
        x_times_2[x] = gfmul(x, 2);
    }

    uint8_t c1 = t1 ^ t4;
    uint8_t c2 = t2 ^ t3 ^ t4;
    uint8_t c3 = t1 ^ t2;
    uint8_t c4 = t1 ^ t3 ^ t4;
    uint8_t c5 = t2 ^ t3;
    uint8_t c6 = t1 ^ t2 ^ t4;
    uint8_t c7 = t3 ^ t4;
    uint8_t c8 = t1 ^ t2 ^ t3;

    for (int case_num = 0; case_num < 4; case_num++) {
        int found = 0;
        results[case_num] = 0;

        for (int i = 0; i < 256; i++) {
            int s = 0;
            for (int x = 0; x < 256; x++) {
                uint8_t condition = 0;
                switch (case_num) {
                    case 0:
                        condition = dot_product(i, x) ^ dot_product(c1, x_times_2[x]) ^ dot_product(c2, x);
                        break;
                    case 1:
                        condition = dot_product(i, x) ^ dot_product(c3, x_times_2[x]) ^ dot_product(c4, x);
                        break;
                    case 2:
                        condition = dot_product(i, x) ^ dot_product(c5, x_times_2[x]) ^ dot_product(c6, x);
                        break;
                    case 3:
                        condition = dot_product(i, x) ^ dot_product(c7, x_times_2[x]) ^ dot_product(c8, x);
                        break;
                }
                if (condition == 0) s++;
                else                s--;
            }
            if (s == 256) {   // perfect correlation found
                results[case_num] = i;
                found = 1;
                break;
            }
        }
    }
}

/* -------------- Hexadecimal string conversion ---------------------- */

/**
 * Convert a hexadecimal string (with or without "0x" prefix) to uint32_t.
 */
uint32_t hex_string_to_uint(const char* hex_str) {
    const char* str = hex_str;
    if (strncmp(str, "0x", 2) == 0 || strncmp(str, "0X", 2) == 0)
        str += 2;

    uint32_t result = 0;
    while (*str) {
        char c = *str;
        uint8_t value;
        if (c >= '0' && c <= '9') value = c - '0';
        else if (c >= 'a' && c <= 'f') value = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') value = c - 'A' + 10;
        else break;
        result = (result << 4) | value;
        str++;
    }
    return result;
}

/* ------------------------- AES S-box -------------------------------- */

static const uint8_t S_BOX[256] = {
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

/* ------------ Pre-computed dot product table ------------------------ */

uint8_t DOT_TABLE[256][256];
int dot_table_computed = 0;

void compute_dot_table() {
    if (dot_table_computed) return;
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            uint8_t result = 0;
            for (int k = 0; k < 8; k++) {
                if (((i >> k) & 1) && ((j >> k) & 1)) result ^= 1;
            }
            DOT_TABLE[i][j] = result;
        }
    }
    dot_table_computed = 1;
}

/* ---------- Structures for correlation statistics ------------------- */

/** Pair of a correlation value and the corresponding l mask. */
typedef struct {
    double correlation;   // signed correlation
    uint8_t l;            // mask l (0..255)
} CorrelationPair;

/** One line (candidate tuple) read from the input file. */
typedef struct {
    char elem2[20];              // hex string of mask h
    char elem3[20];              // hex string of mask xi
    char elem_other[20];         // hex string of another field (unused)
    double elem_value;           // base exponent from e1*e3
    uint8_t results1[4];         // per-byte masks for h
    uint8_t results2[4];         // per-byte masks for xi
} LineData;

/* ----------- Sorting helper for correlations ------------------------ */

int compare_correlation(const void* a, const void* b) {
    const CorrelationPair* pa = (const CorrelationPair*)a;
    const CorrelationPair* pb = (const CorrelationPair*)b;
    double abs_a = fabs(pa->correlation);
    double abs_b = fabs(pb->correlation);
    if (abs_b > abs_a) return 1;
    if (abs_b < abs_a) return -1;
    return 0;
}

/**
 * For a given keystream byte s and a pair (v1, v2) = (h_byte, xi_byte),
 * compute the T best masks l (in terms of |cor|).  The results are stored
 * in s_correlations[t] and s_l_values[t] for t = 0..T-1.
 *
 * The correlation for a given l is computed over all x, y:
 *   cor = (2*valid - total)/total,
 * where valid counts how often
 *   (v1*S(x)) XOR (l*y) XOR (v2*S(z)) = 0
 * with z = s XOR (x+y) mod 256.
 */
void compute_best_l_for_s(uint8_t v1, uint8_t v2, uint8_t s, int T,
                          uint8_t* table_v1, uint8_t* table_v2,
                          double* s_correlations, uint8_t* s_l_values) {
    CorrelationPair all_correlations[256];
    for (int l = 0; l < 256; l++) {
        uint64_t total = 0;
        uint64_t valid = 0;
        for (int x = 0; x < 256; x++) {
            uint8_t v1_sbox = table_v1[x];
            for (int y = 0; y < 256; y++) {
                uint8_t z = s ^ ((x + y) & 0xFF);
                uint8_t v2_sbox = table_v2[z];
                total++;
                if ((v1_sbox ^ DOT_TABLE[l][y] ^ v2_sbox) == 0)
                    valid++;
            }
        }
        double correlation = (2.0 * valid - total) / total;
        all_correlations[l].correlation = correlation;
        all_correlations[l].l = l;
    }
    qsort(all_correlations, 256, sizeof(CorrelationPair), compare_correlation);
    for (int t = 0; t < T; t++) {
        s_correlations[t] = all_correlations[t].correlation;
        s_l_values[t] = all_correlations[t].l;
    }
}

/**
 * For a byte pair (v1, v2) (active bytes of h and xi), compute the
 * T best exponents (log2 of RMS of max |cor|) and also return the
 * (256 x T) matrix of chosen l values (one per s and rank).
 *
 * On exit, *all_zero is set to 1 if all T exponents are essentially zero,
 * and the returned array (of length T) must be freed by the caller.
 */
double* calculate_T_best_exponents_single(uint8_t v1, uint8_t v2, int T,
                                          int* all_zero,
                                          uint8_t*** l_matrix) {
    if (T <= 0 || T > 256) T = 8;
    compute_dot_table();

    uint8_t table_v1[256], table_v2[256];
    for (int i = 0; i < 256; i++) {
        table_v1[i] = DOT_TABLE[v1][S_BOX[i]];
        table_v2[i] = DOT_TABLE[v2][S_BOX[i]];
    }

    double** s_correlations = (double**)malloc(256 * sizeof(double*));
    uint8_t** s_l_values = (uint8_t**)malloc(256 * sizeof(uint8_t*));
    for (int s = 0; s < 256; s++) {
        s_correlations[s] = (double*)malloc(T * sizeof(double));
        s_l_values[s] = (uint8_t*)malloc(T * sizeof(uint8_t));
    }

    for (int s = 0; s < 256; s++) {
        compute_best_l_for_s(v1, v2, s, T, table_v1, table_v2,
                             s_correlations[s], s_l_values[s]);
    }

    double* best_exponents = (double*)malloc(T * sizeof(double));
    *all_zero = 1;
    for (int rank = 0; rank < T; rank++) {
        double sum_abs_correlation_sq = 0.0;
        for (int s = 0; s < 256; s++) {
            double abs_corr = fabs(s_correlations[s][rank]);
            sum_abs_correlation_sq += abs_corr * abs_corr;
        }
        double rms_abs_correlation = sqrt(sum_abs_correlation_sq / 256.0);
        if (rms_abs_correlation > 1e-12) {
            best_exponents[rank] = log2(rms_abs_correlation);
            *all_zero = 0;
        } else {
            best_exponents[rank] = 0.0;
        }
    }

    // second pass: if all exponents are zero, set flag accordingly
    int all_zero_flag = 1;
    for (int i = 0; i < T; i++)
        if (fabs(best_exponents[i]) > 1e-12) { all_zero_flag = 0; break; }
    *all_zero = all_zero_flag;

    *l_matrix = s_l_values;  // caller will free this
    // free the s_correlations storage (not needed outside)
    for (int s = 0; s < 256; s++) free(s_correlations[s]);
    free(s_correlations);
    return best_exponents;
}

/* -------- Recursive generation of combined exponents ---------------- */

/**
 * Generate all combinations, track the maximum sum and the chosen
 * rank (0..7) for each valid byte pair.
 */
void generate_combinations_with_tracking(double elem_value,
                                         double** pair_exponents,
                                         int* valid_indices, int valid_count,
                                         int depth, double current_sum,
                                         int* current_indices,
                                         double* best_sum, int* best_indices) {
    if (depth == valid_count) {
        double sum = elem_value + current_sum;
        if (sum > *best_sum) {
            *best_sum = sum;
            for (int i = 0; i < valid_count; i++)
                best_indices[i] = current_indices[i];
        }
        return;
    }
    int pair_idx = valid_indices[depth];
    for (int t = 0; t < 8; t++) {
        current_indices[depth] = t;
        double new_sum = current_sum + pair_exponents[pair_idx][t];
        generate_combinations_with_tracking(elem_value, pair_exponents,
                                            valid_indices, valid_count,
                                            depth + 1, new_sum,
                                            current_indices, best_sum,
                                            best_indices);
    }
}

/* -------------------- Main processing logic ------------------------- */

/**
 * Read the candidate file, compute conditional correlation exponents,
 * and print the globally best mask combination together with the chosen
 * l bytes.
 */
void process_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }

    char line[1024];
    int line_num = 0;
    LineData* lines = NULL;
    int line_count = 0;

    // Read all lines
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        line[strcspn(line, "\r\n")] = '\0';
        if (strlen(line) == 0) continue;

        char* parts[20];
        int part_count = 0;
        char line_copy[1024];
        strcpy(line_copy, line);
        char* token = strtok(line_copy, " \t");
        while (token != NULL && part_count < 20) {
            parts[part_count++] = token;
            token = strtok(NULL, " \t");
        }
        if (part_count < 5) continue;  // need at least 5 fields

        lines = (LineData*)realloc(lines, (line_count + 1) * sizeof(LineData));
        LineData* ld = &lines[line_count];
        strcpy(ld->elem2, parts[2]);    // mask h
        strcpy(ld->elem3, parts[3]);    // mask xi
        strcpy(ld->elem_other, parts[6]);
        ld->elem_value = atof(parts[6]); // base exponent

        uint32_t val1 = hex_string_to_uint(ld->elem2);
        uint32_t val2 = hex_string_to_uint(ld->elem3);
        allbox_optimized(val1, ld->results1); // per-byte masks for h
        allbox_optimized(val2, ld->results2); // per-byte masks for xi
        line_count++;
    }
    fclose(fp);

    printf("Processing %d lines...\n", line_count);

    double global_best_sum = -1e100;
    int global_best_ranks[4];
    int global_valid_bytes[4];
    uint8_t global_l_values[4];
    int best_line_idx = -1;

    // Process all lines
    for (int idx = 0; idx < line_count; idx++) {
        LineData* ld = &lines[idx];
        double** pair_exponents = (double**)malloc(4 * sizeof(double*));
        uint8_t*** pair_l_matrices = (uint8_t***)malloc(4 * sizeof(uint8_t**));
        int valid_bytes[4];
        int valid_count = 0;

        for (int p = 0; p < 4; p++) {
            uint8_t v1_byte = ld->results1[p];
            uint8_t v2_byte = ld->results2[p];
            int all_zero = 0;
            pair_exponents[p] = calculate_T_best_exponents_single(
                    v1_byte, v2_byte, 8, &all_zero, &pair_l_matrices[p]);
            valid_bytes[p] = all_zero ? 0 : 1;
            if (valid_bytes[p]) valid_count++;
        }

        int* valid_indices = (int*)malloc(valid_count * sizeof(int));
        int idxv = 0;
        for (int i = 0; i < 4; i++) if (valid_bytes[i]) valid_indices[idxv++] = i;

        double best_sum_this_line = -1e100;
        int* best_indices = NULL;
        if (valid_count > 0) {
            best_indices = (int*)malloc(valid_count * sizeof(int));
            int* current_indices = (int*)malloc(valid_count * sizeof(int));
            generate_combinations_with_tracking(ld->elem_value, pair_exponents,
                                                valid_indices, valid_count,
                                                0, 0.0, current_indices,
                                                &best_sum_this_line, best_indices);
            free(current_indices);
        } else {
            best_sum_this_line = ld->elem_value; // no valid bytes, only base
        }

        if (best_sum_this_line > global_best_sum) {
            global_best_sum = best_sum_this_line;
            best_line_idx = idx;
            for (int i = 0; i < 4; i++) global_valid_bytes[i] = valid_bytes[i];
            if (valid_count > 0) {
                for (int i = 0; i < valid_count; i++) {
                    int byte_idx = valid_indices[i];
                    global_best_ranks[byte_idx] = best_indices[i];
                    global_l_values[byte_idx] = pair_l_matrices[byte_idx][0][best_indices[i]];
                }
            }
        }

        if (best_indices) free(best_indices);
        free(valid_indices);
        for (int p = 0; p < 4; p++) {
            free(pair_exponents[p]);
            if (pair_l_matrices[p]) {
                for (int s = 0; s < 256; s++) free(pair_l_matrices[p][s]);
                free(pair_l_matrices[p]);
            }
        }
        free(pair_exponents);
        free(pair_l_matrices);
    }

    // Output the best result
    printf("\n=== Global Maximum Value Details ===\n");
    printf("Global maximum sum: %.6f\n", global_best_sum);
    if (best_line_idx >= 0 && best_line_idx < line_count) {
        LineData* best_line = &lines[best_line_idx];
        printf("Source line index: %d (original line number: %d)\n",
               best_line_idx, best_line_idx + 1);
        printf("Original line data:\n");
        printf("  h   : %s\n", best_line->elem2);
        printf("  xi  : %s\n", best_line->elem3);
        printf("  base exponent : %.6f\n", best_line->elem_value);
        printf("  Per-byte masks (h_byte, xi_byte):\n");
        for (int p = 0; p < 4; p++)
            printf("    Byte %d: h=0x%02x, xi=0x%02x\n",
                   p, best_line->results1[p], best_line->results2[p]);
    }
    printf("\nChosen l values (using s=0 as representative):\n");
    for (int p = 0; p < 4; p++) {
        if (global_valid_bytes[p])
            printf("  Byte %d: l = 0x%02x (rank %d)\n",
                   p, global_l_values[p], global_best_ranks[p]);
        else
            printf("  Byte %d: excluded (all-zero exponents)\n", p);
    }

    free(lines);
}

/* ---------------------------- main ----------------------------------- */

int main() {
    printf("Computation of conditional correlation for SNOW-Vi / SNOW-V\n");
    printf("Only the globally best mask combination is reported.\n");
    printf("================================================================\n\n");

    process_file("0x81ec5a80.txt");

    printf("\nProgram completed successfully.\n");
#ifdef _WIe32
    printf("\nPress any key to exit...");
    getchar();
#endif
    return 0;
}