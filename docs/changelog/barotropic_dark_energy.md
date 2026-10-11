# Barotropic (constant adiabatic sound speed) dark energy model

## Summary

New dark energy model `DarkEnergyBarotropic` (Fortran `TDarkEnergyBarotropic` in
`fortran/DarkEnergyBarotropic.f90`, ini name `barotropic`) for a scalar field with generalized linear
barotropic equation of state `P = ca2*rho + C`, i.e. constant adiabatic sound speed `ca2`, following
[arXiv:1008.1943](https://arxiv.org/abs/1008.1943) (quintessence) and
[arXiv:1206.5194](https://arxiv.org/abs/1206.5194) (phantom). The parameters are `w` (= w0, the value today),
`ca2` and the rest-frame sound speed `cs2` (the `cs2_lam` ini key), giving

```
w(a)            = (1+ca2)(1+w0)/(1 + w0 - (w0 - ca2) a^(3(1+ca2))) - 1
rho(a)/rho(1)   = ((1+w0) a^(-3(1+ca2)) + ca2 - w0)/(1+ca2)
dw/dlna         = 3(1+w)(w - ca2)
```

This is not the CPL `w, wa` parameterization (which CAMB already supports via the `fluid` and `ppf` models):
`ca2` is the asymptotic value of `w` (at early times for `-1 < w0 <= 0`, in the future for phantom `w0 < -1`).

## Implementation notes

- The type extends `TDarkEnergyModel` directly (rather than `TDarkEnergyEqnOfState`) so that the `ca2` field
  sits at a fixed offset visible from Python; the fluid perturbation equations are duplicated from
  `TDarkEnergyFluid` with the exact `dw/dlna` term.
- `w_de` and `grho_de` are written in terms of `(exp(x)-1)/x` (computed accurately for small `x`), so results are
  smooth through `ca2 = -1`, where the closed form has a removable `0/0` and reduces to
  `1 + w = (1+w0)/(1 - 3(1+w0) ln a)`, `rho/rho0 = 1 - 3(1+w0) ln a`.
- Validation: `w0 = -1` is treated as a cosmological constant for any `ca2`; otherwise `w0 <= 0` and `ca2 <= 0` are
  required for `w0 > -1` (else `w > 0` today or at early times), and `ca2 <= w0` for `w0 < -1` (else the density
  changes sign at some `a < 1`). `w` never crosses -1 for valid parameters, so the fluid equations apply.
- `Effective_w_wa` returns `w0` and `-dw/da|_{a=1} = -3(1+w0)(w0-ca2)` for the halofit/HMcode fitting formulae.
- Tested in `testDarkEnergyBarotropic` against the tabulated-`w(a)` fluid model with the same `w(a)`, and against
  the cosmological constant and constant-`w` limits.
