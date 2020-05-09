
const logfact = [0,[sum(log(i) for i=1:n) for n=1:50]...];

"""
    f_logfact(n)

The function returns logarithm of n!
"""
function f_logfact(n)
    (n < 0 || n > 50) && error("n < 0 || n > 50. Modify if needed.")
    @inbounds return logfact[n+1]
end


# special function used in the code:
# (two_n/2)! for even numbers
# (two_n-1/2)! for odd numbers
const logfact2 = [f_logfact(div(two_n,2)) for two_n=1:100];
function f_logfact2(two_n)
    (two_n < 0 || two_n > 100) && error("two_n < 0 || two_n > 100. Modify if needed.")
    @inbounds return logfact2[two_n+1]
end
