# Mask Search for SNOW 5G

This repository contains the codes of the mask search for SNOW 5G presented in Section 4 (Applications to the SNOW-V Family) of our paper, the details being given in the appendix "Supplementary Material for the Search of Linear Masks in SNOW 5G". The search proceeds in four steps, and every step produces the files that are used by the next step. The resulting masks give the highest conditional correlation 2^-62.545 of SNOW 5G.

The codes are listed as follows:

1. Evaluation of the upper bound on cor(P1) for the candidate output masks: the program enumerates the active-byte space and records the correlation of every candidate with its log2 scale and sign (`Correlation computation for the e3.cpp`)
2. Exhaustive evaluation of cor(e1) for all byte pairs (t1, t2), written to `c(e1).txt` (`Correlation computation for the e1.cpp`)
3. Screening of the upper-block masks (n^H, xi^H, h^H) for the fixed mask gamma^H = 0x81ec, the survivors being written to `0x81ec.txt` (`Candidates of (n^H, xi^H, h^H) by cor(e3^H).cpp`)
4. Screening of the lower-block masks (n^L, xi^L, h^L, beta, m) for the fixed mask gamma^L = 0x5a80, the survivors being written to `0x5a80.txt`; the pre-computed cor(e1) is read from `c(e1).txt` (`Candidates of (n^L, xi^L, h^L, beta, m) by cor(e1)cor(e3^L).cpp`)
5. Search of the best byte pair (alpha, l) for every surviving candidate, read from standard input (`Classical_candidates of (alpha, l).cpp`), and evaluation of the classical and of the conditional correlation of the final mask combination (`Computation of classical correlation.cpp`, `Computation of conditional correlation.cpp`)

The programs of steps 3 and 4 use OpenMP. For example, the search is run as follows:

```bash
gcc -O2 -o e3                "Correlation computation for the e3.cpp" -lm
gcc -O2 -o e1                "Correlation computation for the e1.cpp" -lm
gcc -O2 -fopenmp -o screen_H "Candidates of (n^H, xi^H, h^H) by cor(e3^H).cpp" -lm
gcc -O2 -fopenmp -o screen_L "Candidates of (n^L, xi^L, h^L, beta, m) by cor(e1)cor(e3^L).cpp" -lm
gcc -O2 -o classical         "Classical_candidates of (alpha, l).cpp" -lm
gcc -O2 -o clcorr            "Computation of classical correlation.cpp" -lm
gcc -O2 -o cocorr            "Computation of conditional correlation.cpp" -lm

./e3 > "Candidates of gamma by cor(P1).txt"     # step 1
./e1                                            # step 2, asks for a 32-bit hexadecimal mask
./screen_H                                      # step 3, writes 0x81ec.txt
./screen_L                                      # step 4, writes 0x5a80.txt
./classical < 0x81ec.txt                        # step 5
./classical < 0x5a80.txt
./clcorr
./cocorr
```

The file `Candidates of gamma by cor(P1).xlsx` contains the table of step 1 sorted by cor(P1).
