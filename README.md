# 🔬 Quantitative interpretation of lateral adsorbate–intermediate interactions via Butler–Volmer kinetics

**Yilin Zhao, Hengshuo Huang, Mingchuan Luo\***

School of Materials Science and Engineering, Peking University, Beijing 100871, China

*Correspondence*: m.luo@pku.edu.cn

---

## 📌 Overview

This repository contains the MATLAB code used for fitting the modified Butler–Volmer (B–V) kinetic model described in our manuscript. 

The script performs non-linear least-squares fitting of experimental ORR polarization data, extracting three physically meaningful parameters:

| Parameter | Symbol | Unit | Description |
|-----------|--------|------|-------------|
| Exchange current density | `j0` | mA cm⁻² | Intrinsic kinetic rate constant |
| OH energy coefficient | `ε_OH` | V | Descriptor of *OH–intermediate interaction |
| Adsorption equilibrium constant | `K0CA` |

## 📐 Theory

The central kinetic expression (Eq. 13 in the manuscript) for a single adsorbate (*OH) reads:

```
j_K(E) = −j₀ · [1 − θ_OH(E)] · exp(2.303·η / b*) · exp{ε_OH·θ_OH(E) − ln(1 + λ·exp(ε_OH·θ_OH(E)))}
```

where:
- `η = E − E⁰` — overpotential (V vs RHE)
- `b*` — intrinsic Tafel slope (V/dec)
- `θ_OH(E)` — OH adsorption coverage, obtained by cubic spline interpolation of experimental data
- `λ = K0CA` — encapsulates the adsorption equilibrium constant
- `ε_OH` — quantifies the lateral interaction between adsorbed *OH and the reaction intermediate of the rate-determining step


The script also accounts for mass-transport limitations by incorporating the Koutecký–Levich relationship:

```
j_measured = j_K · j_L / (j_K + j_L)
```

where `j_L` is the limiting diffusion current density.


## ⚙️ Requirements

- **MATLAB** R2016a or later (tested on R2023a)
- **Optimization Toolbox** (for `lsqcurvefit`)
- **Statistics and Machine Learning Toolbox** (commented out by default)

## 📊 Data

The script reads experimental data from an Excel file `origin Data-OH.xlsx` containing two sheets:

| Sheet name | Column 1 | Column 2 | Description |
|------------|----------|----------|-------------|
| `sita OH (E)` | Potential (V vs RHE) | θ_OH (monolayer) | OH adsorption coverage vs potential |
| `j vs E` | Current density (mA cm⁻²) | Potential (V vs RHE) | ORR polarization curve |

> ⚠️ **Note:** Place `origin Data-OH.xlsx` in the same directory as the script before running.

## 🚀 Usage

1. Clone or download this repository.
2. Ensure `origin Data-OH.xlsx` is in the working directory.
3. Open and run `Code.m` in MATLAB.

## 📈 Output

The script produces:

1. **Console output** — fitted parameters, convergence status, and residual sum of squares.
2. **Figure 1** — OH coverage vs potential: experimental data points overlaid with the cubic spline interpolation.
3. **Figure 2** — ORR polarization curve: experimental data vs model fit.


