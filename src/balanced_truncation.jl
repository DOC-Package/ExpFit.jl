
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

function baltru_cmplx(a::Vector{ComplexF64}, c::Vector{ComplexF64}, eps::Float64)
    N = length(nu)
    
    # 1. Construct Am = zeros(ComplexF64, N, N) and fill its diagonal with nu.
    Am = zeros(ComplexF64, N, N)
    for i in 1:N
        Am[i,i] = nu[i]
    end

    # 2. Form vectors a and b.
    a_vec = [sqrt(c[i]) for i in 1:N]           # a(i) = c(i)**0.5
    b_vec = [sqrt(conj(c[i])) for i in 1:N]       # b(i) = (conjg(c(i)))**0.5

    # (The Fortran code then deallocates c.)
    # 3. Define x = copy(nu) and y = conj.(nu).
    x = copy(nu)
    y = [conj(nu[i]) for i in 1:N]
    # (nu and c are no longer needed; we’ll replace them later.)

    # 4. Call coneig from the pd_cauchy module.
    #    Assume coneig_cmplx returns a tuple (lam, U).
    lam, U = coneig(a_vec, b_vec, x, y, epsilon)
    # Let n_val be the number of rows of U.
    n_val = size(U, 1)
    

    # (The Fortran block that computes M based on lam is commented out.)
    # We assume here that M (the effective order) is chosen as L.
    M = L

    # 7. Allocate Um as the first M columns of U.
    Um = U[:, 1:M]
    # 8. Compute Umh = herconjg(Um); here we assume this means element‐wise conjugation.
    Umh = conj.(Um)
    
    # 9. Compute Ap = Umh * (Am * conj.(Um))
    Ap = Umh * (Am * conj.(Um))
    # (Deallocate Am and Um by not using them further.)
    
    # 10. Compute bp = Umh * a_vec.
    bp = Umh * a_vec
    # (Deallocate b_vec and Umh.)
    
    # 11. Compute the eigen–decomposition of Ap.
    eig_res = eigen(Ap)
    nu_new = eig_res.values   # new nu (eigenvalues)
    U2 = eig_res.vectors      # eigenvector matrix (M×M)
    
    # 12. Normalize each column of U2 using the sum-of-squares (without conjugation, as in Fortran’s transpose).
    for i in 1:M
        xi = copy(U2[:, i])
        ti = sum(xi .^ 2)
        U2[:, i] .*= 1 / sqrt(ti)
    end

    # 13. Update bp by bp = transpose(U2) * bp.
    bp = transpose(U2) * bp
    # 14. Set new c as bp.^2.
    c_new = [bp[i]^2 for i in 1:M]
    
    
    return a_new, c_new
end
