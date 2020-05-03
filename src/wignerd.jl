
#  _|              _|
#      _|_|_|    _|_|_|_|    _|_|      _|_|_|    _|_|    _|  _|_|
#  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|_|_|  _|_|
#  _|  _|    _|    _|      _|        _|    _|  _|        _|
#  _|  _|    _|      _|_|    _|_|_|    _|_|_|    _|_|_|  _|
#                                          _|
#                                      _|_|

function jacobi_pols(n, a, b, z)
    if (n+a > length(logfact) || n+b > length(logfact))
        error("Error: j is too high, please check the implementation of jacobi polynomials!")
    end

    # special case when I can not calculate log
    if (z ≈ 1) || (z ≈ -1)
        return sum((s % 2 == 0 ? 1.0 : -1.0) *
                exp(logfact[n+a+1] + logfact[n+b+1]-logfact[n-s+1]-logfact[a+s+1]-logfact[s+1]-logfact[n+b-s+1])*
                ((1-z)/2.0)^s*((1+z)/2.0)^(n-s) for s = 0:n)
    end

    # general case
    ls = log((1.0-z)/2.0);
    lc = log((1.0+z)/2.0);
    res = 0.0;
    for s = 0:n
        logs = logfact[(n+a+1)::Int] + logfact[(n+b+1)::Int]-
                    logfact[(n-s+1)::Int]-logfact[(a+s+1)::Int]-
                    logfact[(s+1)::Int]-logfact[(n+b-s+1)::Int];
        args = s*ls + (n-s)*lc;
        res += (s % 2 == 0 ? 1.0 : -1.0) * exp(logs+args);
    end
    return res;
end

function wignerd_hat(j, m1, m2, z)
        if abs(m1) > j || abs(m2) > j
            return zero(z)
        end
        factor = ((abs(m1-m2)+m1-m2)/2) % 2 == 0 ? one(z) : -one(z);
        am1 = abs(m1); am2 = abs(m2);
        M = (am1 > am2) ? am1 : am2;
        N = (am1 < am2) ? am1 : am2;
        gammas = logfact[Int(j-M+1)]+logfact[Int(j+M+1)]-(logfact[Int(j-N+1)]+logfact[Int(j+N+1)])
        return factor / 2^M * exp(gammas/2)*
               jacobi_pols(Int(j-M), Int(abs(m1-m2)), Int(abs(m1+m2)), z);
end

"""
    wignerd(j, m1, m2, cosβ)

Small wigner d-function for representation `j` with indices `m1`, `m2`, the argument `cosβ` is the cosine of the rotation angle.
The function gives the value of the matrix element
```
    ⟨ j m1 | exp(-i Jy β)| j m2 ⟩.
```
The input values are expected to be __integers__.
For a general case including half-integers see `wignerd_doublearg`.
"""
function wignerd(j, m1, m2, z)
    (z ≈ 1) && return m1 == m2 ? one(z) : zero(z)
    (z ≈ -1) && return m1 == -m2 ? (iseven(j-m2) ? one(z) : -one(z)) : zero(z)
    #
    hat = wignerd_hat(j, m1, m2, z);
    xi = sqrt(1-z)^abs(Int(m1-m2))*sqrt(1+z)^abs(Int(m1+m2));
    return hat*xi;
end

"""
    wignerD(j, m1, m2, α, cosβ, γ)

Wigner D-function for representation `j` with indices `m1`, `m2`, α, β, and γ are the rotation angles
The function gives the value of the matrix element
```
    ⟨ j m1 | exp(-i Jz α) exp(-i Jy β) exp(-i Jz γ) | j m2 ⟩.
```
The input values are expected to be __integers__.
For a general case including half-integers see `wignerD_doublearg`.
"""
wignerD(j, m1, m2, α, cosβ, γ)  = wignerd(j, m1, m2, cosβ) * cis(-m1*α-m2*γ)

#
#  _|                  _|      _|_|              _|              _|
#  _|_|_|      _|_|_|  _|    _|                      _|_|_|    _|_|_|_|    _|_|      _|_|_|    _|_|    _|  _|_|
#  _|    _|  _|    _|  _|  _|_|_|_|  _|_|_|_|_|  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|_|_|  _|_|
#  _|    _|  _|    _|  _|    _|                  _|  _|    _|    _|      _|        _|    _|  _|        _|
#  _|    _|    _|_|_|  _|    _|                  _|  _|    _|      _|_|    _|_|_|    _|_|_|    _|_|_|  _|
#                                                                                        _|
#                                                                                    _|_|

function wignerd_hat_doublearg(two_j, two_m1, two_m2, z)
        # @show j, m1, m2
        if abs(two_m1) > two_j || abs(two_m2) > two_j
            return zero(z)
        end
        factor = (abs(two_m1-two_m2)+two_m1-two_m2) % 8 == 4 ? -one(z) : one(z);
        two_am1 = abs(two_m1); two_am2 = abs(two_m2);
        (two_M, two_N) = (two_am1 > two_am2) ? (two_am1,two_am2) : (two_am2,two_am1);
        #
        j_mnus_M = div(two_j-two_M,2)
        j_plus_M = div(two_j+two_M,2)
        j_mnus_N = div(two_j-two_N,2)
        j_plus_N = div(two_j+two_N,2)
        #
        gammas = logfact[j_mnus_M+1] +
                 logfact[j_plus_M+1] -
                 logfact[j_mnus_N+1] -
                 logfact[j_plus_N+1]
        return factor / 2^(two_M/2) * exp(gammas/2) *
               jacobi_pols(j_mnus_M, div(abs(two_m1-two_m2),2), div(abs(two_m1+two_m2),2), z);
end

"""
    wignerd_doublearg(j, m1, m2, cosβ)

Small wigner d-function for representation `j` with indices `m1`, `m2`, the argument `cosβ` is the cosine of the rotation angle.
The function gives the value of the matrix element
```
    ⟨ j m1 | exp(-i Jy β)| j m2 ⟩.
```
The function requires __doubled value of the momenta__ for the input.
"""
function wignerd_doublearg(two_j, two_m1, two_m2, z)
    (z ≈ 1) && return two_m1 == two_m2 ? one(z) : zero(z)
    (z ≈ -1) && return two_m1 == -two_m2 ? (iseven(div(two_j-two_m2,2)) ? one(z) : -one(z)) : zero(z)
    #
    hat = wignerd_hat_doublearg(two_j, two_m1, two_m2, z);
    xi = (1-z)^(abs(two_m1-two_m2)/4)*(1+z)^(abs(two_m1+two_m2)/4);
    return hat*xi;
end

"""
    wignerD_doublearg(two_j, two_m1, two_m2, α, cosβ, γ)

Wigner D-function for representation `j` with indices `m1`, `m2`, α, β, and γ are the rotation angles
The function gives the value of the matrix element
```
    ⟨ j m1 | exp(-i Jz α) exp(-i Jy β) exp(-i Jz γ) | j m2 ⟩.
```
The function requires __doubled value of the momenta__ for the input.
"""
wignerD_doublearg(two_j, two_m1, two_m2, α, cosβ, γ) =
    wignerd_doublearg(two_j, two_m1, two_m2, cosβ) * cis(-two_m1*α/2-two_m2*γ/2);

"""
    kronecker(i, j) = (i==j) ? 1 : 0
"""
kronecker(i, j) = (i==j) ? 1 : 0
