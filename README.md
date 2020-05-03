# PartialWaveFunctions

[![Build Status](https://travis-ci.com/mmikhasenko/PartialWaveFunctions.jl.svg?branch=master)](https://travis-ci.com/mmikhasenko/PartialWaveFunctions.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mmikhasenko/PartialWaveFunctions.jl?svg=true)](https://ci.appveyor.com/project/mmikhasenko/PartialWaveFunctions-jl)
[![Codecov](https://codecov.io/gh/mmikhasenko/PartialWaveFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/PartialWaveFunctions.jl)

Julia native implementation of the special functions used in the Partial Wave Analysis for High Energy Physics. Currently, Wigner D-functions and Clebsch-Gordan coefficients are available.

## Installation
```julia
] add https://github.com/mmikhasenko/PartialWaveFunctions.jl # to be registered soon
```

## Usage
```julia
using PartialWaveFunctions

# convenient call for integer indices
let j=3, (m1,m2) = (1,-1), cosθ=0.3
   wignerd(j,m1,m2,cosθ)
end # return 0.293

clebschgordan(1,0,1,0,1,0) #  # <1, 0; 1, 0 | 1, 0> = 0.0 : ρ⁰ → π⁰ π⁰
CG(1,0,1,0,1,0) # a shortcut
```

general implementation including half-integer indices
```julia
#
let two_j=3, (two_m1,two_m2) = (1,-1), cosθ=0.3
   wignerd_doublearg(two_j,two_m1,two_m2, cosθ)
end # return -0.562

clebschgordan_doublearg(2,0,1,1,1,1) # <1, 0; 1/2, 1/2 | 1/2, 1/2> = -0.577
CG_doublearg(2,0,1,1,1,1) # a shortcut

```

## Similar packages:
 * python calls via `SymPy.jl`. Ideal for symbolic calculations. Works pretty with jupyter notebooks due to the latex output. See details in the [test/physics](https://github.com/JuliaPy/SymPy.jl/blob/master/test/test-physics.jl).
 * [WignerD.jl](https://github.com/jishnub/WignerD.jl) interface `Fortran` for the `WignerD`.
 * [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl) `Julia` package specialized on Symbols. Particularly it contains the Clebsch-Gordan symbols.
 * [GSL.jl](https://github.com/JuliaMath/GSL.jl) interface `C++`. It can calculate Sperical Harmionics, Legendre polynomials. `WignerD` is not [wrapper-up](https://github.com/JuliaMath/GSL.jl/issues/66).
