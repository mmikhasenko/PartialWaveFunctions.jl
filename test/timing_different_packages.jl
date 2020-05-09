using Test
using Parameters
using BenchmarkTools

#
import PartialWaveFunctions: CG_doublearg, clebschgordan_doublearg
ClGd_pwf(two_j1,two_m1,two_j2,two_m2,two_j,two_m) = clebschgordan_doublearg(two_j1,two_m1,two_j2,two_m2,two_j,two_m)

#
using GSL
function ClGd_gsl(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    factor = sqrt(two_j+1)*(mod(two_j1-two_j2+two_m,4)==2 ? -1 : +1)
    three_j = sf_coupling_3j(two_j1,two_j2,two_j,two_m1,two_m2,-two_m)
    return factor*three_j;
end

#
using SymPy
import PyCall
PyCall.pyimport_conda("sympy.physics.wigner","sympy")
import_from(sympy.physics.wigner)
#
ClGd_sympy(two_j1,two_m1,two_j2,two_m2,two_j,two_m) =
    convert(Float64, clebsch_gordan(Sym(two_j1)/2, Sym(two_j2)/2, Sym(two_j)/2, Sym(two_m1)/2, Sym(two_m2)/2, Sym(two_m)/2))

#
import HalfIntegers: half
using WignerSymbols
function ClGd_WS(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    ((abs(two_m1) > two_j1) || (abs(two_m2) > two_j2) || (abs(two_m ) > two_j )) && return 0.0
    return clebschgordan(half(two_j1),half(two_m1),half(two_j2),half(two_m2),half(two_j),half(two_m))
end

#                                      _|
#  _|  _|_|    _|_|_|  _|_|_|      _|_|_|
#  _|_|      _|    _|  _|    _|  _|    _|
#  _|        _|    _|  _|    _|  _|    _|
#  _|          _|_|_|  _|    _|    _|_|_|


function rand_clebsch(;two_j_max::Int=15)
    two_j1 = rand(0:two_j_max); two_m1 = rand(-two_j1:2:two_j1)
    two_j2 = rand(0:two_j_max); two_m2 = rand(-two_j2:2:two_j2)
    two_j = rand(abs(two_j2-two_j1):2:(two_j1+two_j2))
    two_m = two_m1+two_m2
    abs(two_m) > two_j && return rand_clebsch(; two_j_max=two_j_max)
    (two_j1=two_j1, two_j2=two_j2, two_j=two_j,
     two_m1=two_m1, two_m2=two_m2, two_m=two_m)
end

# fraction of vanishing clebshces ≈ 0.6%
sum(let
    @unpack two_j1, two_j2, two_j, two_m1, two_m2, two_m = rand_clebsch()
    ClGd_pwf(two_j1,two_m1,two_j2,two_m2,two_j,two_m) ≈ 0
end for _ in 1:10_000)


#                                                                    _|
#    _|_|_|    _|_|    _|_|_|  _|_|    _|_|_|      _|_|_|  _|  _|_|        _|_|_|    _|_|    _|_|_|
#  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  _|_|      _|  _|_|      _|    _|  _|    _|
#  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  _|        _|      _|_|  _|    _|  _|    _|
#    _|_|_|    _|_|    _|    _|    _|  _|_|_|      _|_|_|  _|        _|  _|_|_|      _|_|    _|    _|
#                                      _|
#                                      _|

for _ in 1:1_000
    @unpack two_j1, two_j2, two_j,
            two_m1, two_m2, two_m = rand_clebsch()
    #
    v1 = ClGd_pwf(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    v2 = ClGd_gsl(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    v3 = ClGd_sympy(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    v4 = ClGd_WS(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    #
    (abs(v1-v2) > 1e-5) && error("gsl is different: $(v1) $(v2)")
    (abs(v1-v3) > 1e-5) && error("sympy is different")
    (abs(v1-v4) > 1e-5) && error("WS is different")
end

#    _|      _|                  _|
#  _|_|_|_|      _|_|_|  _|_|        _|_|_|      _|_|_|
#    _|      _|  _|    _|    _|  _|  _|    _|  _|    _|
#    _|      _|  _|    _|    _|  _|  _|    _|  _|    _|
#      _|_|  _|  _|    _|    _|  _|  _|    _|    _|_|_|
#                                                    _|
#                                                _|_|

@btime for _ in 1:1_000
    @unpack two_j1, two_j2, two_j, two_m1, two_m2, two_m = rand_clebsch()
    ClGd_pwf(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
end # 227.200 μs

@btime for _ in 1:1_000
    @unpack two_j1, two_j2, two_j, two_m1, two_m2, two_m = rand_clebsch()
    ClGd_gsl(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
end # 921.8 μs

@btime for _ in 1:1_000
    @unpack two_j1, two_j2, two_j, two_m1, two_m2, two_m = rand_clebsch()
    ClGd_WS(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
end # 7.483 ms

@btime for _ in 1:1_000
    @unpack two_j1, two_j2, two_j, two_m1, two_m2, two_m = rand_clebsch()
    ClGd_sympy(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
end # 1.61 s
