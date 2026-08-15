# On Skew Cyclic Codes Over $\mathbb{Z}_4 + u\mathbb{Z}_4$: Structure, Gray Images, and New Quaternary Linear Codes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

Official Python implementation and computational search engine accompanying the research manuscript:
> **"On Skew Cyclic Codes Over $\mathbb{Z}_4 + u\mathbb{Z}_4$: Structure, Gray Images, and New Quaternary Linear Codes"**  
> *Basri Çalışkan, Erhan Turan, and Cennet Eskal*

---

## 📌 Overview

Exhaustive search over polynomial rings is computationally intractable due to exponential state space expansion ($\mathcal{O}(16^n)$). This framework implements a **two-phase algebraically-guided search engine** that collapses the effective combinatorial search domain to $\mathcal{O}(16^{n/2})$:

1. **Symbolic Phase (SymPy):** Enforces quotient relations ($u^2 = 2$), computes the ring automorphism $\sigma$, and strictly restricts candidate generation to self-dual structural invariants ($c_i = -\theta(c_{n-1-i})$).
2. **Vectorized Numerical Funnel:** Maps candidate structures via vectorized Gray maps, performs Gaussian elimination over $\mathbb{Z}_4$, and filters candidates via rank-profile pruning prior to high-cost minimum Lee distance ($d_L$) evaluations.
3. **Automated Magma Cross-Validation:** Features a pure-Python client that verifies candidate weight spectra and target Lee distances via the online Magma calculator API—requiring **no local Magma installation**.

---

## 🚀 Key Results & Benchmarks

| Length $n$ | Gray Length $2n$ | Full Search Space | Constrained Space | Best $d_L$ Found | Avg. Iterations / Runtime |
|:---:|:---:|:---:|:---:|:---:|:---:|
| **$n=8$** | $16$ | $\approx 4.29 \times 10^9$ | $65,536$ | **$8$** | $\approx 100\text{ iter } (\sim 83\text{ s})$ |
| **$n=8$** | $16$ | $\approx 4.29 \times 10^9$ | $65,536$ | **$8,4$** | $\approx 100\text{ iter } (\sim 93\text{ s})$ |
| **$n=16$** | $32$ | $\approx 1.84 \times 10^{19}$ | $\approx 4.29 \times 10^9$ | **$8$** | $\approx 10,000\text{ iter } (\sim 6\text{ h } 05\text{ m})$ |

*Benchmarks executed on an Intel® Core™ i7-13620H processor with 16 GB RAM.*

---

## 🛠️ Installation & Setup

```bash
# Clone the repository
git clone [https://github.com/drerhanturan/codegenerator.git](https://github.com/drerhanturan/codegenerator.git)
cd codegenerator

# Install dependencies
pip install -r requirements.txt
```

## 🚀 Usage

The framework is organized modularly via the CodeGenerator class, allowing easy integration into custom search pipelines or verification scripts.

1. Stochastic Guided Search
Configure structural ring maps and execute heuristic search runs:

```python
from CodeGenerator import CodeGenerator

# Initialize search engine
generator = CodeGenerator()

# Define automorphism theta, derivations, and Gray mapping components
thetas = ['a+2*b', 'b']
derivations = ['2*b', '2*b']
grays = ['3*b', 'a+b']

# Set parameters: u constraint, polynomial degree / length parameter x, base ring Z_4
generator.set_parameters(u=1, x=7, z=4, thetas=thetas, derivations=derivations, grays=grays)

# Run guided search: (target_found_count=10, max_iterations=100)
# Halts when 10 new/better codes are discovered or iteration threshold is reached
generator.solve_for_random_polynomials(10, 100)

```
2. Single Polynomial Evaluation
Evaluate and debug specific candidate polynomials directly:

```python
from CodeGenerator import CodeGenerator

generator = CodeGenerator()
generator.set_parameters(u=1, x=7, z=4, thetas=['a+2*b', 'b'], derivations=['2*b', '2*b'], grays=['3*b', 'a+b'])

# Define polynomial coefficients: e.g., 3x^5 + (3u+1)x^4 + x^3 + 3ux^2 + (3u+2)x + (3u+1)
coeffs = ['3', '3*u+1', '1', '3*u', '3*u+2', '3*u+1']

result = generator.solve_for_custom_polynomial(coeffs)
print(result.get_debug_info())
print("Linearity Condition:", result.condition)
```

## 📁 Output Artifacts
Benchmark runs dynamically benchmark discovered codes against known tables (quantumcodes.info $\mathbb{Z}_4$ database) and partition results into two summary text files:

- usefuls.txt: Contains novel (no previous code known for the given $[n, 4^{k_1} 2^{k_2}]$ profile) and better codes (strictly exceeding previously known minimum Lee distances $d_L$).
- equals.txt: Contains discovered codes that match existing best-known minimum Lee distance bounds.



## 📖 Citation

This repository provides the official implementation and reproducibility assets for our research manuscript. If you use the algorithms, ring constructions, or benchmark results in your work, please cite the underlying paper:

```bibtex
@unpublished{caliskan2026skewz4,
  author = {Basri \c{C}al{\i}\c{s}kan and Erhan Turan and Cennet Eskal},
  title  = {On Skew Cyclic Codes Over $\mathbb{Z}_4 + u\mathbb{Z}_4$: Structure, Gray Images, and New Quaternary Linear Codes},
  note   = {Manuscript submitted for publication},
  year   = {2026}
}
```
