#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/*
 * This program evaluates the correlation of a specific linear mask combination
 * for the SNOW 5G stream cipher.  It uses the modular addition automaton to
 * compute sub-correlations and then combines them according to the structure
 * of the third sub-approximation (e3) in the theoretical linear cryptanalysis.
 *
 * The main function enumerates the two active bytes of the intermediate mask theta
 * and sums the products of two pairs of automaton
 * outputs:
 *   Upper block (e3^H):
 *     automaton(theta^H, n^H, gamma^H) * automaton(xi^H, h^H, theta)
 *     with n^H = 0x81e8, gamma^H = 0x81ec, xi^H = 0x8180, h^H = 0xe1a0.
 *
 *   Lower block (e3^L):
 *     automaton(theta^L, n^L, gamma^L) * automaton((xi \oplus beta)^L, h^L, theta)
 *     with n^L = 0x5a00, gamma^L = 0x5a80, (xi \oplus beta) ^L = 0x7f00, h^L = 0x7f40.
 *
 * The correlated product is then printed together with the sign.
 * All correlations are printed in log2 scale (absolute value).
 */

/* Automaton for 16-bit modular addition 
 *  input1, input2 : input masks
 *  output         : output mask
 *  Returns the correlation (signed) as a double, or -100.0 for invalid state
*/
double automaton16(uint32_t input1, uint32_t input2, uint32_t output)
{
    double k = 0;          // exponent counter (solid arrows)
    int e = 0;             // automaton state (0 or 1)
    int t;
    int flag;
    double sign = 1.0;     // sign due to transitions that flip correlation

    for (t = 0; t < 16; t++)
    {
        // Form the transition label mu_i = u_i*4 + v_i*2 + w_i
        flag = ((output >> (15 - t)) & 1) * 4
             + ((input1 >> (15 - t)) & 1) * 2
             + ((input2 >> (15 - t)) & 1);

        if (e == 0)
        {
            if (flag == 7)
            {
                e = 1;
            }
            else if (flag == 0)
            {
                e = 0;
            }
            else
            {
                return -100.0;   // non-zero correlation impossible
            }
        }
        else   // e == 1
        {
            if ((flag == 0) || (flag == 3) || (flag == 5) || (flag == 6))
            {
                e = 1;
                k++;             // solid arrow while staying in e1
            }
            else
            {
                e = 0;
                k++;             // solid arrow while moving from e1 to e0
            }
            if ((flag == 3) || (flag == 4))
            {
                sign = (-1.0) * sign;   // sign change
            }
        }
    }
    return sign * pow(2.0, (-1) * k);
}

int main()
{
    uint32_t theta;
    double cor_h = 0.0;   // correlation sum for the upper block
    double cor_l = 0.0;   // correlation sum for the lower block
    double cor;
    int i;

    // Enumerate theta = (i << 8) for i = 0..255
    // Upper block (e3^H)
    for (i = 0; i < 256; i++)
    {
        theta = i << 8;
        double prod_h1 = automaton16(theta, 0x81e8, 0x81ec);
        double prod_h2 = automaton16(0x8180, 0xe1a0, theta);
        if ((prod_h1 > -100) && (prod_h2 > -100))
            cor_h += prod_h1 * prod_h2;
    }
    printf("Upper block cor (log2|cor|): %.4f\n", log2(fabs(cor_h)));

    // Lower block (e3^L)
    for (i = 0; i < 256; i++)
    {
        theta = i << 8;
        double prod_l1 = automaton16(theta, 0x5a00, 0x5a80);
        double prod_l2 = automaton16(0x7f00, 0x7f40, theta);
        if ((prod_l1 > -100) && (prod_l2 > -100))
            cor_l += prod_l1 * prod_l2;
    }
    printf("Lower block cor (log2|cor|): %.4f\n", log2(fabs(cor_l)));

    cor = cor_h * cor_l;
    printf("Total cor (log2|cor|): %.5f\n", log2(fabs(cor)));
    if (cor > 0.0)
        printf("Sign: positive (+1)\n");
    else
        printf("Sign: negative (-1)\n");

    return 0;
}