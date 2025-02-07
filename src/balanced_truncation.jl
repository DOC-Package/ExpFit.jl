
#
# 1. Build Am = diag(nu) ∈ ℂ^(N×N)
# 2. Form a vector a by taking element‐wise square roots of c and a vector b
#    by taking the square roots of the conjugates of c.
# 3. Let x = nu and y = conj.(nu); then “deallocate” the original nu and c.
# 4. Call coneig (from the pd_cauchy module) with (N, a, b, x, y, epsilon) to get
#    lam and U.
# 5. For each column of U, perform a phase–adjustment:
#       ti = sum(U(:,i).^2)
#       phase = ti / abs(ti)
#       scale = sqrt(conj(phase))
#       U(:,i) *= scale
# 6. (A block that would choose M based on lam is commented out in Fortran.)
#    Instead, we assume M is given (or chosen as the input L parameter).
# 7. Form Um as the first M columns of U.
# 8. Compute Umh = conj.(Um) [Fortran’s herconjg(Um)]
# 9. Compute Ap = Umh * (Am * conj.(Um))
# 10. Compute bp = Umh * a.
# 11. Compute the eigen–decomposition of Ap, so that
#       eigen(Ap) returns eigenvalues (to be stored in nu) and eigenvectors U.
# 12. For each column i of the eigenvector matrix U, normalize it as:
#       ti = sum(U[:,i].^2)
#       scale = 1 / sqrt(ti)
#       U[:,i] *= scale
# 13. Replace bp by (transpose(U) * bp) and then set
#       c[i] = (bp[i])^2 for i = 1:M.
# 14. Finally, sort (using a tandem sort routine) alphaout and rhoout.

function balanced_truncation(a::Vector{ComplexF64}, c::Vector{ComplexF64}, eps::Float64)
    N = length(nu)

    # 2. Form vectors a and b.
    #a_vec = [sqrt(c[i]) for i in 1:N]           # a(i) = c(i)**0.5
    #b_vec = [sqrt(conj(c[i])) for i in 1:N]       # b(i) = (conjg(c(i)))**0.5

    # (The Fortran code then deallocates c.)
    # 3. Define x = copy(nu) and y = conj.(nu).
    #x = copy(nu)
    #y = [conj(nu[i]) for i in 1:N]

    sv, U = coneig(sqrt.(c), sqrt.(conj(c)), nu, conj(nu))
    M = findfirst(i -> 2 * sum(sv[i+1:end]) < eps, 1:N)

    # 7. Allocate Um as the first M columns of U.
    Um = U[:, 1:M]
    
    #Am = zeros(ComplexF64, N, N)
    #for i in 1:N
    #    Am[i,i] = nu[i]
    #end
    # 9. Compute Ap = Umh * (Am * conj.(Um))
    Ap = Um' * (Diagonal(nu) * conj.(Um))
    
    bp = Umh * sqrt.(c)
    
    # 11. Compute the eigen–decomposition of Ap.
    res = eigen(Ap)
    a_new = res.values  
    V = res.vectors
    
    # 12. Normalize each column of U2 using the sum-of-squares (without conjugation, as in Fortran’s transpose).
    for i in 1:M
        xi = copy(V[:, i])
        ti = sum(xi .^ 2)
        V[:, i] .*= 1 / sqrt(ti)
    end
    bp = transpose(V) * bp
    c_new = bp.^2
    
    return a_new, c_new
end
