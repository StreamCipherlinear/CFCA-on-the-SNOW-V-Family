# Conditional Fast Correlation Attacks with Applications to the SNOW-V Family

This repository contains the source code of our paper "Conditional Fast Correlation Attacks with Applications to the SNOW-V Family". It implements the search of the conditional linear masks and the computation of the correlations (Section 4, Applications to the SNOW-V Family), and the validation of the framework on a reduced variant of SNOW 5G (Section 5, Experimental Verification).

The codes are organized as follows:

* The mask search and the correlation computation for SNOW 5G and for SNOW-V / SNOW-Vi are placed in the directory `(Section 4) Applications to the SNOW-V Family`
* The validation of the framework on the reduced variant of SNOW 5G is placed in the directory `(Section 5) Experimental Verification`

Each directory contains its own ReadMe that describes the programs and how to run them. All programs are plain C/C++ and need no external library; GCC or Clang is recommended, and the programs that use OpenMP are compiled with `-fopenmp`.

## Licence

This code is released under the MIT Licence, see `LICENSE`. If you use it in your work, please cite the corresponding paper.
