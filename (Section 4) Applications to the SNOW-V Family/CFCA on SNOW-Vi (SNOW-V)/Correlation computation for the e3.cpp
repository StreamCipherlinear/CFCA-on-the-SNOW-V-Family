/*
 * Evaluation of the correlation of e3 for SNOW-Vi (and SNOW-V).
 *
 * This program computes the correlation contributed by the third sub-approximation
 * e3 for a specific linear mask combination.  It uses the 32-bit modular addition
 * automaton (+ mod 2^32) to evaluate the contributions of the two sub-approximations
 * P1 and P2 and then sums their products over an enumerated intermediate mask theta.
 *
 * The fixed masks used in this example are:
 *   n    = 0x81ec5a00
 *   gamma    = 0x81ec5a80
 *   xi    = 0x81800000
 *   h    = 0x81c00000
 *
 * The intermediate mask theta is enumerated as theta = i << 8 for i = 0..255.
 * This corresponds to the active byte of theta being in the second byte position
 * (the remaining bytes are zero), consistent with the linear mask template
 * derived for SNOW-Vi.
 *
 * The correlation is computed as:
 *   cor(e3) = Sum_theta automaton(theta, n, gamma) * automaton(xi, h, theta)
 *
 * The output prints log2(|cor|) and the sign (+1 or -1).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/**
 * 32-bit modular addition automaton (+ mod 2^32).
 *
 * @param input1   first input mask  (v)
 * @param input2   second input mask (w)
 * @param output   output mask       (u)
 * @return signed correlation value, or -100.0 if zero correlation is forced.
 */
double automaton32(uint32_t input1, uint32_t input2, uint32_t output)
{
    double k = 0;          // exponent counter (solid-arrow transitions)
    int e = 0;             // automaton state (0 or 1)
    int t;
    int flag;
    double sign = 1.0;     // sign factor

    for (t = 0; t < 32; t++)
    {
        // Form the transition label mu_i = u_i*4 + v_i*2 + w_i
        flag = ((output >> (31 - t)) & 1) * 4
             + ((input1 >> (31 - t)) & 1) * 2
             + ((input2 >> (31 - t)) & 1);

        if (e == 0)
        {
            if (flag == 7)
                e = 1;
            else if (flag == 0)
                e = 0;
            else
                return -100.0;   // non-zero correlation impossible
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
                sign = -sign;    // sign change
        }
    }
    return sign * pow(2.0, -k);
}

int main()
{
    uint32_t theta;
    double cor = 0.0;           // correlation accumulator (must be initialised)
    int i;

    // Enumerate the active byte of theta (second byte, i = 0..255)
    for (i = 0; i < 256; i++)
    {
        theta = i << 8;         // active byte in position 1, others zero

        double c1 = automaton32(theta, 0x81ec5a00, 0x81ec5a80);   // cor(P1)
        double c2 = automaton32(0x81800000, 0x81c00000, theta);   // cor(P2)

        if ((c1 > -100.0) && (c2 > -100.0))
            cor += c1 * c2;
    }

    printf("log2(|cor|): %.5f\n", log2(fabs(cor)));
    if (cor > 0.0)
        printf("Sign: +1\n");
    else
        printf("Sign: -1\n");

    return 0;
}