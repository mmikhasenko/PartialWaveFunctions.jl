# PartialWaveFunctions

[![Build Status](https://travis-ci.com/mmikhasenko/PartialWaveFunctions.jl.svg?branch=master)](https://travis-ci.com/mmikhasenko/PartialWaveFunctions.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mmikhasenko/PartialWaveFunctions.jl?svg=true)](https://ci.appveyor.com/project/mmikhasenko/PartialWaveFunctions-jl)
[![Codecov](https://codecov.io/gh/mmikhasenko/PartialWaveFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/PartialWaveFunctions.jl)

Julia native implementation of the special functions used in the Partial Wave Analysis for High Energy Physics. Currently, Wigner D-functions and Clebsch-Gordan coefficients are available.

## Installation
```julia
] add PartialWaveFunctions
```

## Usage
```julia
using PartialWaveFunctions

# convenient call for integer indices
let j=3, (m1,m2) = (1,-1), cosθ=0.3
   wignerd(j,m1,m2,cosθ)
end # return 0.293

clebschgordan(1,0,1,0,1,0) # <1, 0; 1, 0 | 1, 0> = 0.0 : ρ⁰ → π⁰ π⁰
CG(1,0,1,0,1,0) # a shortcut
```

General implementation includes the half-integer indices:
```julia
let two_j=3, (two_m1,two_m2) = (1,-1), cosθ=0.3
   wignerd_doublearg(two_j,two_m1,two_m2, cosθ)
end # return -0.562

clebschgordan_doublearg(2,0,1,1,1,1) # <1, 0; 1/2, 1/2 | 1/2, 1/2> = -0.577
CG_doublearg(2,0,1,1,1,1) # a shortcut
```

## Related packages:
 * python calls via `SymPy.jl`. Ideal for symbolic calculations. Works pretty with jupyter notebooks due to the latex output. See details in the [test/physics](https://github.com/JuliaPy/SymPy.jl/blob/master/test/test-physics.jl).
 * [WignerD.jl](https://github.com/jishnub/WignerD.jl) interfaces `Fortran` for the `WignerD`.
 * [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl) `Julia` package specialized on Symbols. Particularly it contains the Clebsch-Gordan coefficients.
 * [GSL.jl](https://github.com/JuliaMath/GSL.jl) interfaces `C++`. It can calculate Sperical Harmionics, Legendre polynomials. `WignerD` is not [wrapped-up](https://github.com/JuliaMath/GSL.jl/issues/66).

## References
 * The Wigner functions are expressed via the Jacobi polynomials Pₙ⁽ᵃᵇ⁾(z) using Eq. (3.74) of
    L. Biedenharn, J. Louck, and P. Carruthers, Angular Momentum in Quantum Physics: Theory and Application
 * The Jacobi polynomials Pₙ⁽ᵃᵇ⁾(z) are codded using a series expression in powers of (1-z), see e.g. [wikipedia page](https://en.wikipedia.org/wiki/Jacobi_polynomials).
 * Clebsch-Gordan coefficients are computed from explicit expression via a finite series, see e.g. [wikipedia page](https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients)
