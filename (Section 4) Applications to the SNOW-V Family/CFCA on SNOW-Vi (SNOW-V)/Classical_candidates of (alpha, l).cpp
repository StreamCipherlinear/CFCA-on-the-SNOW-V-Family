/*
 *
 * Search for optimal (alpha, l) for each (h, xi) byte pair in the classical
 * correlation evaluation of e2 -- for SNOW-V and SNOW-Vi.
 *
 * This program reads pre-screened mask candidates from standard input
 * (typically redirected from the file "0x81ec5a80.txt" produced by the
 * e3 filtering step) and, for every active-byte position of each tuple,
 * exhaustively searches all non-zero byte pairs (alpha, l) to maximise
 * the contribution of e2.  The overall classical correlation exponent
 * is then obtained by summing the best per-byte exponents together with
 * the base exponent from the previous e1*e3 stage.  The line (i.e. mask
 * combination) with the highest total exponent is reported.
 *
 * Input format (one line per tuple, 8 fields):
 *   table_idx  mask_n  mask_h  mask_xi  mask_beta  mask_m  exponent_base  unused
 * The seventh field (0-based index 6) is the base exponent (log2|cor|)
 * coming from the e1*e3 product.  Fields 2 and 3 are the 32-bit masks
 * h and xi (hexadecimal, with or without "0x" prefix).
 *
 * The masks h and xi are first decomposed into their four active bytes
 * via the allbox_optimized() function.  For each of these four byte
 * positions, all (alpha, l) in {1..255}^2 are tested, and the pair that
 * gives the highest |correlation| is retained.  The program finally
 * outputs the optimal (alpha, l) for each byte position of the best
 * overall candidate.
 *
 * Usage:
 *   ./Classical_candidates_alpha_l_SNOW_V < 0x81ec5a80.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* ---------------------------- GF(2^8) helpers ---------------------------- */

/**
 * Multiplication in GF(2^8) with the AES polynomial x^8+x^4+x^3+x+1.
 */
uint8_t gfmul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    uint8_t hi_bit_set;
    for (int i = 0; i < 8; i++) {
        if (b & 1) p ^= a;
        hi_bit_set = a & 0x80;
        a <<= 1;
        if (hi_bit_set) a ^= 0x1b;
        b >>= 1;
    }
    return p;
}

/**
 * Dot product (= parity of bitwise AND), returns 0 or 1.
 */
uint8_t dot_product(uint8_t a, uint8_t b) {
    uint8_t result = a & b;
    result ^= (result >> 4);
    result ^= (result >> 2);
    result ^= (result >> 1);
    return result & 0x01;
}

/**
 * For a 32-bit mask t = (t4,t3,t2,t1) (big-endian byte order),
 * find the 8-bit output masks that give perfect correlation (+/-256)
 * for the four AES S-box equations (cases 0..3), if they exist.
 * Results are stored in results[0..3]; a zero value indicates none found.
 */
void allbox_optimized(uint32_t t, uint8_t* results) {
    uint8_t t1 = t & 0xFF;
    uint8_t t2 = (t >> 8) & 0xFF;
    uint8_t t3 = (t >> 16) & 0xFF;
    uint8_t t4 = (t >> 24) & 0xFF;

    uint8_t x_times_2[256];
    for (int x = 0; x < 256; x++) x_times_2[x] = gfmul(x, 2);

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
                    case 0: condition = dot_product(i, x) ^ dot_product(c1, x_times_2[x]) ^ dot_product(c2, x); break;
                    case 1: condition = dot_product(i, x) ^ dot_product(c3, x_times_2[x]) ^ dot_product(c4, x); break;
                    case 2: condition = dot_product(i, x) ^ dot_product(c5, x_times_2[x]) ^ dot_product(c6, x); break;
                    case 3: condition = dot_product(i, x) ^ dot_product(c7, x_times_2[x]) ^ dot_product(c8, x); break;
                }
                if (condition == 0) s++; else s--;
            }
            if (s == 256) { results[case_num] = i; found = 1; break; }
        }
    }
}

/* ------------------ Hexadecimal string to integer ------------------- */

uint32_t hex_string_to_uint(const char* hex_str) {
    const char* str = hex_str;
    if (strncmp(str, "0x", 2) == 0 || strncmp(str, "0X", 2) == 0) str += 2;
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

/* --------------------------- AES S-box ------------------------------- */

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

/* --------- Pre-computed dot product table (parity of bitwise AND) --- */

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

/* --- Classical correlation evaluation for one (h, xi) byte pair ----- */

/**
 * For fixed (h, xi) at a single active-byte position, this function
 * evaluates the e2 contribution for a given (alpha, l) and returns
 * log2(|cor|).  If the correlation is zero, *all_zero is set to 1.
 *
 * The computation follows the exact e2 formula:
 *   cor = (1/2^24) * Sum_{x,y,z} (-1)^{xi*S(y) XOR alpha*(y+z) XOR l*z XOR alpha*x XOR h*S(x)}.
 */
double compute_exponent_for_pair(uint8_t h, uint8_t xi, uint8_t alpha, uint8_t l, int* all_zero) {
    compute_dot_table();

    // Step 1: Count for all (y,z) the value of T = xi*S(y) XOR alpha*(y+z) XOR l*z
    uint64_t count0 = 0, count1 = 0;
    for (int y = 0; y < 256; y++) {
        uint8_t xi_sy = DOT_TABLE[xi][S_BOX[y]];
        for (int z = 0; z < 256; z++) {
            uint8_t alpha_y_plus_z = DOT_TABLE[alpha][(y + z) & 0xFF];
            uint8_t l_z = DOT_TABLE[l][z];
            uint8_t T = xi_sy ^ alpha_y_plus_z ^ l_z;
            if (T == 0) count0++; else count1++;
        }
    }

    // Step 2: Count total valid = sum over x of (U=0? count0 : count1)
    uint64_t valid = 0;
    for (int x = 0; x < 256; x++) {
        uint8_t alpha_x = DOT_TABLE[alpha][x];
        uint8_t h_sx = DOT_TABLE[h][S_BOX[x]];
        uint8_t U = alpha_x ^ h_sx;
        if (U == 0) valid += count0;
        else        valid += count1;
    }

    uint64_t total = 256LL * 256 * 256;  // 16777216
    double corr = (2.0 * valid - total) / total;
    double abs_corr = fabs(corr);

    if (abs_corr > 1e-12) {
        *all_zero = 0;
        return log2(abs_corr);
    } else {
        *all_zero = 1;
        return 0.0;
    }
}

/* -------------- Data structure for one input line ------------------- */

typedef struct {
    double elem_value;          // base exponent from e1*e3
    uint8_t results1[4];        // per-byte masks for h
    uint8_t results2[4];        // per-byte masks for xi
    char elem2_str[20];         // hex string for mask h (debug)
    char elem3_str[20];         // hex string for mask xi (debug)
    char original_line[1024];   // copy of the whole input line
} LineData;

/* ----------------------- Main processing logic ---------------------- */

void process_file() {
    LineData* lines = NULL;
    int line_count = 0;
    char line[1024];

    printf("Please enter data lines (format: ... elem2 elem3 ... elem_value ...), Ctrl+D to end:\n");
    while (fgets(line, sizeof(line), stdin)) {
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

        lines = (LineData*)realloc(lines, (line_count + 1) * sizeof(LineData));
        LineData* ld = &lines[line_count];

        // Save original line for later display
        strncpy(ld->original_line, line, sizeof(ld->original_line) - 1);
        ld->original_line[sizeof(ld->original_line) - 1] = '\0';

        // Parse: fields 2 and 3 are h and xi, field 6 is the base exponent
        strcpy(ld->elem2_str, parts[2]);
        strcpy(ld->elem3_str, parts[3]);
        ld->elem_value = atof(parts[6]);

        uint32_t val1 = hex_string_to_uint(ld->elem2_str);
        uint32_t val2 = hex_string_to_uint(ld->elem3_str);

        allbox_optimized(val1, ld->results1); // h bytes
        allbox_optimized(val2, ld->results2); // xi bytes

        line_count++;
    }

    if (line_count == 0) {
        printf("No valid data read.\n");
        free(lines);
        return;
    }

    printf("Read %d lines of data.\n", line_count);

    // Total iterations for progress reporting
    int total_pairs = line_count * 4;
    int total_alpha_l = 255 * 255;
    long long total_iterations = (long long)total_pairs * total_alpha_l;
    long long iter_done = 0;
    int progress_step = total_iterations / 100;
    if (progress_step == 0) progress_step = 1;

    // Allocate per-pair best (alpha, l) storage
    double (*pair_best_exp)[4] = (double(*)[4])malloc(line_count * 4 * sizeof(double));
    int (*pair_best_alpha)[4] = (int(*)[4])malloc(line_count * 4 * sizeof(int));
    int (*pair_best_l)[4] = (int(*)[4])malloc(line_count * 4 * sizeof(int));
    int (*pair_has_valid)[4] = (int(*)[4])calloc(line_count, 4 * sizeof(int));

    if (!pair_best_exp || !pair_best_alpha || !pair_best_l || !pair_has_valid) {
        printf("Memory allocation failed!\n");
        free(lines);
        free(pair_best_exp);
        free(pair_best_alpha);
        free(pair_best_l);
        free(pair_has_valid);
        return;
    }

    // Initialise best exponents to negative infinity
    for (int i = 0; i < line_count; i++)
        for (int p = 0; p < 4; p++)
            pair_best_exp[i][p] = -1e100;

    printf("Searching best (alpha,l) for each byte pair...\n");

    for (int row = 0; row < line_count; row++) {
        for (int p = 0; p < 4; p++) {
            uint8_t h_byte = lines[row].results1[p];
            uint8_t xi_byte = lines[row].results2[p];

            double best_exp = -1e100;
            int best_alpha = 0, best_l = 0, has_valid = 0;

            for (int alpha = 1; alpha < 256; alpha++) {
                for (int l = 1; l < 256; l++) {
                    int all_zero;
                    double exp_val = compute_exponent_for_pair(h_byte, xi_byte, alpha, l, &all_zero);
                    iter_done++;
                    if (iter_done % progress_step == 0) {
                        printf("Progress: %.2f%%\r", 100.0 * iter_done / total_iterations);
                        fflush(stdout);
                    }

                    if (!all_zero && exp_val > best_exp) {
                        best_exp = exp_val;
                        best_alpha = alpha;
                        best_l = l;
                        has_valid = 1;
                    }
                }
            }

            pair_best_exp[row][p]   = best_exp;
            pair_best_alpha[row][p] = best_alpha;
            pair_best_l[row][p]     = best_l;
            pair_has_valid[row][p]  = has_valid;
        }
    }

    printf("\n\n"); // finish progress line

    // Determine the row with the highest total exponent
    double best_standard_value = -1e100;
    int best_row = -1;
    for (int row = 0; row < line_count; row++) {
        double sum_exp = 0.0;
        for (int p = 0; p < 4; p++)
            if (pair_has_valid[row][p])
                sum_exp += pair_best_exp[row][p];
        double standard_value = lines[row].elem_value + sum_exp;
        if (standard_value > best_standard_value) {
            best_standard_value = standard_value;
            best_row = row;
        }
    }

    if (best_row != -1) {
        printf("=== Optimal row by standard value ===\n");
        printf("Standard value: %f\n", best_standard_value);
        printf("Original line: %s\n", lines[best_row].original_line);
        printf("Best (alpha, l) for each byte position:\n");
        for (int p = 0; p < 4; p++) {
            if (pair_has_valid[best_row][p]) {
                printf("  Byte %d: alpha = 0x%02x, l = 0x%02x, exp = %f\n",
                       p, pair_best_alpha[best_row][p], pair_best_l[best_row][p],
                       pair_best_exp[best_row][p]);
            } else {
                printf("  Byte %d: no valid exponent found\n", p);
            }
        }
    } else {
        printf("No valid data rows found.\n");
    }

    free(pair_best_exp);
    free(pair_best_alpha);
    free(pair_best_l);
    free(pair_has_valid);
    free(lines);
}

int main() {
    printf("Classical (alpha,l) search for e2 (SNOW-V / SNOW-Vi)\n");
    printf("Data is read from standard input, each line format (8 fields):\n");
    printf("  ... <h> <xi> ... ... ... <base_exp> ...\n");
    printf("Please enter data (Ctrl+D to end).\n\n");
    process_file();
    printf("\nProgram finished.\n");
    return 0;
}