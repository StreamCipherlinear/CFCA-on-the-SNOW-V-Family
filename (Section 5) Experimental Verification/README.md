# Experimental Verification of the Conditional Fast Correlation Attack

This repository contains the codes for the validation presented in Section 5 (Experimental Verification) of our paper. All programs are self-contained C/C++ programs that print their results to stdout and need no input files.

The codes are all based on the reduced variant of SNOW 5G used in the paper: the LFSR is 64 bits wide with the primitive polynomial f(x) = x^64 + x^4 + x^3 + x + 1 and advances 16 steps per clock, the two 16-bit taps T1 and T2 are taken from bit positions 32-47 and 0-15, and the three 128-bit registers R1, R2 and R3 are reduced to 16-bit words, each viewed as a 2x2 matrix of 4-bit nibbles. The FSM uses S_small, the 4-bit S-box of [musa2003simplified] applied to each nibble followed by a MixColumn over GF(16), in place of the AES round function of SNOW 5G, as done in [fastsnow2]. The registers are updated by

```
R1^{t+1} = R2^t + (R3^t ^ T2^t)          (16-bit modular addition)
R2^{t+1} = S_small(R1^t)
R3^{t+1} = S_small(R2^t)
```

and the keystream word is z^t = (R1^t + T1^t) ^ R2^t.

Outer linear masks alpha, beta, gamma are applied to the keystream words z^{t-1}, z^t, z^{t+1}, and linear masks l, m, n, h to the taps T1^{t-1}, T1^t, T1^{t+1}, T2^t. With an intermediate mask xi on R2, the approximation is partitioned into the three sub-approximations e1, e2 and e3, with

```
cor(e) = cor(e1) * Sum_xi cor(e2) * cor(e3).
```

The codes are listed as follows:

* The masks of the sub-approximations e1 and e2, where the program reports |cor(e1)| = 0.25 and |cor(e2)| = 0.125 (`ce1e2.cpp`)
* The optimal combination of (beta, gamma, xi, h, n) maximizing the classical correlation, c = 2^-1.586 (`ce3.cpp`)
* The classical correlation c(N) of the three-word approximation, 2^-6.211, which is the optimal classical correlation reported in Section 5 (`cor_traverse_xi.cpp`)
* The adaptive (l, tau) for each u, and the theoretical conditional correlation with the RMS over u (`cor_con.cpp`)
* The independence between the sub-approximations e2 and e3 in the classical case (`Q1_KL_e2e3_classical.cpp`) and in the conditional case (`Q1_KL_e2e3_conditional.cpp`), through the KL divergence between the joint distribution and the product of the marginals
* The WHT relation between the classical and conditional correlations stated in Proposition 1, which prints True if all checks pass (`Q2_Verify_proposition1.cpp`)
* The accuracy of the conditional correlation estimated from parity-check equations: the theoretical conditional correlation 2^-6.0924, and M = 2^18.53 without collisions and 2^30.24 with collisions (`Q3_Verify_conditional_correlation_two_mode.cpp`)
* The recovery of the initial LFSR state: with 2^18.53 keystream blocks the correct initial state ranks first (`Q4_Recover_initial_state.cpp`)

In the codes, `e1`, `e2`, `e3` denote the sub-noises $e_1$, $e_2$, $e_3$ of the paper, `xi` is the intermediate mask, `tau` the sign mask, `u` the conditioning value, `B` the unknown bits, `K` the collision-key bits, and `M` the number of parity-check equations.

Q3 and Q4 are C++17: Q4 builds with GCC 6 or later, whereas Q3 uses structured bindings and therefore needs GCC 7 or later (or Clang 4 or later); both also build with MSVC (`cl /O2 /std:c++17 /EHsc`). The other codes use the GCC/Clang builtin `__builtin_parity` and need GCC or Clang. Q1_* use OpenMP and are the heaviest runs. For example:

```bash
g++ -O2 -o ce1e2            ce1e2.cpp -lm
g++ -O2 -o ce3               ce3.cpp -lm
g++ -O2 -o cor_con           cor_con.cpp -lm
g++ -O2 -o cor_traverse_xi   cor_traverse_xi.cpp -lm
g++ -O2 -fopenmp -o Q1_classical   Q1_KL_e2e3_classical.cpp -lm
g++ -O2 -fopenmp -o Q1_conditional Q1_KL_e2e3_conditional.cpp -lm
g++ -O2 -o Q2                Q2_Verify_proposition1.cpp -lm
g++ -O2 -std=c++17 -o Q3     Q3_Verify_conditional_correlation_two_mode.cpp -lm
g++ -O2 -std=c++17 -o Q4     Q4_Recover_initial_state.cpp -lm

./ce1e2                     # Optimal m = 0x08, |cor(e1)| = 0.250000
./ce3                        # c = +0.333008 (sign=+1, log2(|c|)=-1.586)
./cor_traverse_xi            # c(N) = +0.013498306, log2(|c(N)|) = -6.211078
./cor_con                    # RMS of |cor| over u = 0.153093, log2(RMS) = -2.707519
./Q1_classical               # Maximum KL divergence: D_KL = ...
./Q2                         # True
./Q3 --symbols 20000 --repetitions 2 --seed 1     # Verdict (3-sigma): PASS
./Q4 --mode direct --repetitions 2 --seed 1       # rank=1/65536, SUCCESS
```

Q3 takes the options `--seed`, `--symbols`, `--pair-symbols`, `--repetitions`, `--B`, `--L`, `--K` and `--ps` (defaults: B = 16, L = 24, K = 16, ps = 0.999992, repetitions = 100), and Q4 takes `--mode direct|pair`, `--B`, `--K`, `--symbols`, `--equations`, `--ps`, `--seed`, `--repetitions` and `--help`.
