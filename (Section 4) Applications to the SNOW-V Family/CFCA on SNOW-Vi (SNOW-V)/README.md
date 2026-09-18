# Mask Search for SNOW-V and SNOW-Vi

This repository contains the codes of the mask search for SNOW-V and SNOW-Vi presented in Section 4 (Applications to the SNOW-V Family) of our paper, the details being given in the appendix "Supplementary Material for the Search of Linear Masks in SNOW-Vi" and "Linear Mask Search of SNOW-Vi". The search proceeds in four steps, and every step produces the files that are used by the next step. The resulting masks give the highest conditional correlation 2^-42.999 of SNOW-V and SNOW-Vi.

The codes are listed as follows:

1. Evaluation of the upper bound on cor(P1) for the candidate output masks: the program enumerates the active-byte space and prints the correlation of every candidate with its log2 scale and sign (`Correlation computation for the e3.cpp`)
2. Exhaustive evaluation of cor(e1) for all byte pairs (t1, t2), written to `c(e1).txt` (`Correlation computation for the e1.cpp`)
3. Screening of the masks (n, xi, h, beta, m) by cor(e1)*cor(e3), where the 32-bit modular addition is handled by its automaton; the pre-computed cor(e1) is read from `c(e1).txt` and the survivors are written to `0x81ec5a80.txt` (`Candidates of (n, xi, h, beta, m) by cor(e1)cor(e3).cpp`)
4. Search of the best byte pair (alpha, l) for every surviving candidate, read from standard input (`Classical_candidates of (alpha, l).cpp`), and evaluation of the conditional correlation of the final mask combination, read from `0x81ec5a80.txt` (`Computation of conditional correlation.cpp`)

The program of step 3 uses OpenMP. For example, the search is run as follows:

```bash
gcc -O2 -o e3        "Correlation computation for the e3.cpp" -lm
gcc -O2 -o e1        "Correlation computation for the e1.cpp" -lm
gcc -O2 -fopenmp -o screen "Candidates of (n, xi, h, beta, m) by cor(e1)cor(e3).cpp" -lm
gcc -O2 -o classical "Classical_candidates of (alpha, l).cpp" -lm
gcc -O2 -o cocorr    "Computation of conditional correlation.cpp" -lm

./e3 > "Candidates of gamma by cor(P1).txt"     # step 1
./e1                                            # step 2, asks for a 32-bit hexadecimal mask
./screen                                        # step 3, reads c(e1).txt, writes 0x81ec5a80.txt
./classical < 0x81ec5a80.txt                    # step 4
./cocorr
```
