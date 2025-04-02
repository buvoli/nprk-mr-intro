---

This directory contains MATLAB code for reproducing the numerical experiment (Figure 4) from the manuscript

    T. Buvoli, B. K. Tran, B. S. Southworth, "Multirate Runge-Kutta for Nonlinearly Partitioned Systems", 2025.

All code is released under the BSD license (See LICENSE).

---

Directory Structure: 

| directory     |  description |
| :--           | :--          |
| figures       | All figures are generated in this directory |
| integrators   | <ol> <li> General implementations of NPRK and ERK methods </li> <li> Functions to generate NPRK coefficient tensors </li> <li> MR-NPRK implementations that follow manuscript pseudocodes. </li> </ol> |
| problems      | <ol> <li> Burgers with nonlinear diffusion and spatial-varying coefficient </li> <li> Viscous Burgers </li> </ol>

**Running Code**: Figure 4 in the manuscript can be reproduced by running script "main_svbnd_paper.m"

**Important:** Set variable `use_cache = false` to generate new timing data.

---
