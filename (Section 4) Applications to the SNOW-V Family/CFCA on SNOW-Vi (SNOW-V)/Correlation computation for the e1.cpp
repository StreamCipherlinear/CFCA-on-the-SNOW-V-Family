/*
 * AES S-box Linear Correlation Analysis for e1 (SNOW-V family)
 *
 * This program computes the correlation of the first sub-approximation e1
 * for the SNOW-V, SNOW-Vi, and SNOW 5G stream ciphers.
 *
 * Given a 32-bit hexadecimal input representing a combination of linear
 * masks, the program evaluates the S-box related part of e1 over all
 * possible values of the associated state bytes.  It uses precomputed tables
 * for the AES S-box, parity, and Galois-field multiplication by 2 in order
 * to accelerate the computation.
 *
 * The results, consisting of the log2-scaled absolute correlation and the
 * sign for each (t1, t2) pair, are printed to the console and saved to
 * the file c(e1).txt.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

// AES S-box
static const uint8_t S[256] = {
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

// Precomputed parity table for all 8-bit values
static uint8_t parity_table[256];

// Precomputed multiplication by 2 in GF(2^8)
static uint8_t mul2_table[256];

/**
 * @brief Initialise lookup tables for performance.
 */
void init_tables(void) {
    // Initialise parity table
    for (int i = 0; i < 256; i++) {
        uint8_t parity = 0;
        uint8_t temp = i;
        for (int j = 0; j < 8; j++) {
            parity ^= (temp & 1);
            temp >>= 1;
        }
        parity_table[i] = parity;
    }
    
    // Initialise multiplication-by-2 table in GF(2^8)
    for (int i = 0; i < 256; i++) {
        uint8_t a = i;
        uint8_t carry = a & 0x80;
        a <<= 1;
        if (carry) a ^= 0x1B;
        mul2_table[i] = a;
    }
}

/**
 * @brief Fast parity computation using the precomputed table.
 */
uint8_t fast_parity(uint8_t value) {
    return parity_table[value];
}

/**
 * @brief Fast multiplication by 2 in GF(2^8) using the precomputed table.
 */
uint8_t fast_mul2(uint8_t value) {
    return mul2_table[value];
}

/**
 * @brief Fast dot product (parity of bitwise AND) using the parity table.
 */
uint8_t fast_dot_product(uint8_t a, uint8_t b) {
    return fast_parity(a & b);
}

/**
 * @brief Optimised correlation analysis for a given 32-bit input.
 *
 * The input encodes four byte-masks that are used to evaluate the
 * S-box related linear approximation of e1.  The correlation values
 * (log2 scale and sign) are printed for every (t1, t2) pair and
 * written to c(e1).txt.
 */
void optimized_correlation_analysis(uint32_t hex_input, FILE *file_out) {
    uint8_t t1_byte = hex_input & 0xFF;
    uint8_t t2_byte = (hex_input >> 8) & 0xFF;
    uint8_t t3_byte = (hex_input >> 16) & 0xFF;
    uint8_t t4_byte = hex_input >> 24;
    
    printf("Input bytes: 0x%02X, 0x%02X, 0x%02X, 0x%02X\n\n", 
           t4_byte, t3_byte, t2_byte, t1_byte);
    printf("t1\tt2\tCorrelation (log2)\tSign\n");
    printf("-----------------------------------------------\n");
    
    // Precompute the S-box linear approximation for every input x
    uint8_t linear_approx[256];
    for (int x = 0; x < 256; x++) {
        uint8_t sbox_out = S[x];
        uint8_t mul2_sbox = fast_mul2(sbox_out);
        linear_approx[x] = 
            fast_dot_product(t1_byte ^ t4_byte, mul2_sbox) ^
            fast_dot_product(t2_byte ^ t3_byte ^ t4_byte, sbox_out);
    }
    
    int num = 0;
    for (int t1 = 0; t1 < 256; t1++) {
        for (int t2 = 0; t2 < 256; t2++) {
            int32_t correlation_sum = 0;
            
            for (int x = 0; x < 256; x++) {
                uint8_t lin_approx_val = linear_approx[x];
                
                for (int y = 0; y < 256; y++) {
                    uint8_t z = (x + y) & 0xFF;
                    uint8_t condition = 
                        fast_dot_product(t1, z) ^ 
                        lin_approx_val ^ 
                        fast_dot_product(t2, y);
                    
                    correlation_sum += (condition == 0) ? 1 : -1;
                }
            }
            
            // Output only significant correlations
            if (abs(correlation_sum) >= 512) {
                double log_corr = log2(abs(correlation_sum));
                int sign = (correlation_sum > 0) ? 1 : -1;
                
                printf("0x%02x\t0x%02x\t%.6f\t\t%d\n", 
                       t1, t2, log_corr - 16, sign);
                fprintf(file_out, "0x%02x 0x%02x %.6f \t%d\n", 
                        t1, t2, log_corr - 16, sign);
                num++;
            }    
        }
    }
    
    printf("num = %d\n", num);
    fprintf(file_out, "num = %d\n", num);
}

int main() {
    uint32_t hex_input;
    
    init_tables();
    
    printf("========== AES S-box Linear Correlation Analysis ==========\n\n");
    printf("Optimised version with precomputed tables\n\n");
    printf("Enter 32 bit hexadecimal input (e.g., 0x12345678): ");
    b
    if (scanf("0x%x", &hex_input) != 1) {
        printf("Invalid input format. Please use format: 0x12345678\n");
        return 1;
    }
    
    printf("\n");
    
    FILE *file_out = fopen("c(e1).txt", "w");
    if (file_out == NULL) {
        printf("Error: Cannot open output file c(e1).txt\n");
        return 1;
    }
    
    clock_t start_time = clock();
    optimized_correlation_analysis(hex_input, file_out);
    clock_t end_time = clock();
    
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    printf("\nExecution time: %.2f seconds\n", execution_time);
    fprintf(file_out, "\nExecution time: %.2f seconds\n", execution_time);
    
    fclose(file_out);
    printf("\nResults have been written to c(e1).txt\n");
    
    return 0;
}