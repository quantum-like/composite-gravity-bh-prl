# Composite-gravity black hole — PRL Letter data and code

Symbolic checks and figure assets for the Physical Review Letters submission on the static composite-gravity black hole (draft title in the scripts: “Gauge-regular but metric-singular”). This repository is the data and code availability package for that submission. It ships the verification programs as provided and the merged invariant figure. It does not re-derive the theory.

Conventions stated in the scripts follow Öttinger, *J. Phys. Commun.* **10**, 105003 (2026): tetrad (45), metric (46), solution (51), composition rule (4), field strength (18), and field equations (19). Units in the scripts are $r_0 = 1$ (geometrized, $r_0 = Gm/c^2$).

## How to run

Python 3 with [SymPy](https://www.sympy.org/) and [mpmath](https://mpmath.org/). NumPy and Matplotlib are not imported; the figure is already rendered.

```bash
pip install -r requirements.txt
python verify_composite_bh_checks.py
python verify_grok_independent.py
python verify_49_solves_19.py
```

`verify_composite_bh_checks.py` is the exact rational SymPy check. `verify_grok_independent.py` does not import it; it repeats the connection, $F\cdot F$, and field-equation comparisons with mpmath at 80 decimal digits and 6th-order differences, and a full run takes several minutes. Those two scripts print to stdout and do not write files. `verify_49_solves_19.py` prints to stdout and writes `results/verify_49_vs_19_out.txt`, creating `results/` if it is missing.

## What the scripts check

### `verify_composite_bh_checks.py`

Exact checks. Cartesian coordinates with $r = |\mathbf{x}|$, evaluated on Pythagorean lattice points so that $r$ is rational, plus a radial reduction $z = \sqrt{r^2 - x^2 - y^2}$ for identities in $r$.

- **Radial reduction.** Coordinate-condition residuals of the tetrad, the condition labeled (32), after that radial substitution, and the on-axis component $A_{(01)0}$ from the composition rule (4).
- **(49) versus (19).** Pointwise residual between the composition rule (4) and the closed form (49), reported as $\max |A_{(4)} - A_{(49)}|$ on the lattice (points with $r = r_0 = 1$ are skipped). The field strength is then built from (49) and the field equation (19) is evaluated two ways: the same commutator sign, and the opposite sign. The script prints the max absolute residual on the lattice and whether the consistent-sign residual is identically zero after the radial reduction.
- **$F \cdot F$.** The Lie-algebra invariant $\sum_a F_a \cdot F^a$ at three lattice radii, for the consistent commutator sign and the opposite sign, next to the draft formula $4(6r^2 + 4r + 1) / [r^4 (r + 1)^4]$.
- **PPN $\beta = \gamma = 1$.** Weak-field series, in $u = r_0/R$, of the metric functions
  $B(R) = \bigl[(2R - 3)/(2R - 1)\bigr]^2$ and
  $g_{RR}(R) = \bigl[(2R + 1)/(2R - 1)\bigr]^2$.
  The printed coefficients are the post-Newtonian comparison ($\beta = \gamma = 1$ is the Schwarzschild leading structure).
- **Photon sphere ($R$) and Killing surface gravity ($\kappa$).** Roots of $\partial_R(R^2/B) = 0$, printed as photon-sphere values of $R/r_0$, and the Killing surface gravity
  $\kappa(R) = \partial_R\sqrt{B} \,/\, \sqrt{g_{RR}}$,
  including the value at $R = 3/2$.

### `verify_grok_independent.py`

Independent of the SymPy script (it does not import it). The connection is the Levi-Civita spin connection of tetrad (45)+(51), in matrix form. Derivatives are 6th-order central differences. Spacetime indices are raised with $\eta = \mathrm{diag}(-1,1,1,1)$.

- **(49).** At six Cartesian points, time and spatial connection entries are printed next to the (49) prefactor $1/[r^2(r+1)]$ times the corresponding coordinate.
- **$F \cdot F$.** $- \tfrac{1}{2} \sum_{m,n} \mathrm{tr}(F_{mn} F^{mn})$ next to the same draft formula as above.
- **(19) residual.** $\max |D_m F^{m\nu}|$ for $\nu = 0,1,2,3$ with the commutator sign used in $F_{mn} = \partial_m w_n - \partial_n w_m + [w_m, w_n]$, and one opposite-sign divergence component ($\nu = 1$) so the two signs can be told apart.

### `verify_49_solves_19.py`

Dedicated (49)-versus-(19) residual check. SymPy builds the static solution (49), the field strength (18), and the Yang–Mills residual of field equation (19) for both commutator signs. mpmath samples that residual at Cartesian points with $r_0 = 1$, including a point with $r = r_0$. The script also prints the Lorenz-gauge residual (38).

## Figures

| File | Contents |
| --- | --- |
| `figures/FIG1_merged_invariants_FF_B.pdf` | Merged invariant figure ($F \cdot F$ and $B$), vector |
| `figures/FIG1_merged_invariants_FF_B.png` | Same figure, raster |

These are the `paper1_figures_rev1` assets supplied with the scripts. The verification scripts do not regenerate them.

## License

[MIT License](LICENSE). Copyright 2026 Amir Vahid, Ocean Reviver Inc.

## Citation

```bibtex
@misc{vahid2026compositegravitybh,
  author       = {Amir Vahid},
  title        = {Symbolic checks for a composite-gravity static black hole},
  year         = {2026},
  organization = {Ocean Reviver Inc.},
  howpublished = {GitHub repository accompanying a Physical Review Letters submission},
  url          = {https://github.com/quantum-like/composite-gravity-bh-prl}
}
```

Amir Vahid, Ocean Reviver Inc. This code and figure set accompanies the PRL submission; it is not a substitute for the Letter.
