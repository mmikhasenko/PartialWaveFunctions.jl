"""
    CG(j1,m1,j2,m2,j,m)

a shortcut for `clebschgordan`
"""
CG(j1,m1,j2,m2,j,m) = clebschgordan(j1,m1,j2,m2,j,m)


"""
    CG_doublearg(j1,m1,j2,m2,j,m)

a shortcut for `clebschgordan_doublearg`
"""
CG_doublearg(two_j1,two_m1,two_j2,two_m2,two_j,two_m) = clebschgordan_doublearg(two_j1,two_m1,two_j2,two_m2,two_j,two_m)


"""
    clebschgordan(j1,m1,j2,m2,j,m)

gives a numerical value of the Clebsch-Gordan coefficient
```
    ⟨ j₁ m₁ ; j₂ m₂ | j m ⟩
```
The input values are expected to be __integers__.
For a general case including half-integers see `clebschgordan_doublearg`.
"""
clebschgordan(j1,m1,j2,m2,j,m) = clebschgordan_doublearg(2j1,2m1,2j2,2m2,2j,2m)


"""
    clebschgordan_doublearg(two_j1,two_m1,two_j2,two_m2,two_j,two_m)

gives a numerical value of the Clebsch-Gordan coefficient
```
    ⟨ j₁ m₁ ; j₂ m₂ | j m ⟩
```
The function requires __doubled value of the momenta__ for the input.
"""
function clebschgordan_doublearg(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    ((abs(two_m1) > two_j1) || (abs(two_m2) > two_j2) || (abs(two_m ) > two_j )) && return 0.0
    ((two_m1+two_m2 != two_m) || !(abs(two_j1-two_j2) ≤ two_j ≤ two_j1+two_j2))  && return 0.0
     #
    prefactor = sqrt(two_j+1)*
        exp( ( f_logfact2(two_j1+two_j2-two_j) +
               f_logfact2(two_j1+two_j -two_j2) +
               f_logfact2(two_j2+two_j -two_j1) -
               f_logfact2(two_j1+two_j2+two_j+2) +
               #
               f_logfact2(two_j1+two_m1) +
               f_logfact2(two_j1-two_m1) +
               f_logfact2(two_j2+two_m2) +
               f_logfact2(two_j2-two_m2) +
               f_logfact2(two_j +two_m ) +
               f_logfact2(two_j -two_m ) ) / 2)
    res = 0.0
    two_t_min = max(0,
                    two_j2-two_m1-two_j,
                    two_j1+two_m2-two_j)
    two_t_max = min(two_j1+two_j2-two_j,
                    two_j1-two_m1,
                    two_j2+two_m2)
    #
    for two_t = two_t_min:2:two_t_max
        logs = f_logfact2(two_t) +
               f_logfact2(two_j-two_j2+two_m1+two_t) +
               f_logfact2(two_j-two_j1-two_m2+two_t) +
               f_logfact2(two_j1+two_j2-two_j-two_t) +
               f_logfact2(two_j1-two_m1-two_t) +
               f_logfact2(two_j2+two_m2-two_t);
        res += (abs(two_t) % 4 == 2 ? -1.0 : 1.0) * exp(-logs);
    end
    res *= prefactor
    return res;
end
