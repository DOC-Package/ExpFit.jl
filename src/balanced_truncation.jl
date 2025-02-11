function balanced_truncation_sub(a::AbstractVector{T}, c::AbstractVector{T}, Um::AbstractMatrix{ComplexF64}) where T<:Number

    Ap = Um' * (Diagonal(a) * conj.(Um))
    bp = Um' * sqrt.(c)
    
    res = eigen(Ap)
    a_new = res.values  
    V = res.vectors
    for i in 1:size(Um, 2)
        xi = copy(V[:, i])
        ti = sum(xi .^ 2)
        V[:, i] .*= 1 / sqrt(ti)
    end
    bp = transpose(V) * bp
    c_new = bp.^2
    
    return Exponentials(a_new, c_new)
end

"""
    balanced_truncation(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, eps::Float64) :: Exponentials

Given the exponents `a` and coefficients `c` of a sum of exponentials, compute a new sum of exponentials with a reduced number of terms for a given tolerance `eps`.
"""
function balanced_truncation(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, eps::Float64)
    a = ComplexF64.(a); c = ComplexF64.(c)
    sv, U = coneig(sqrt.(c), sqrt.(conj.(c)), a, conj.(a))
    M = findfirst(i -> 2 * sum(sv[i+1:end]) < eps, 1:length(sv))
    Um = U[:, 1:M]
    return balanced_truncation_sub(a, c, Um)
end

"""
    balanced_truncation(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, M::Int) :: Exponentials

Given the exponents `a` and coefficients `c` of a sum of exponentials, compute a new sum of exponentials with a given number of terms `M`.
"""
function balanced_truncation(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, M::Int)
    a = ComplexF64.(a); c = ComplexF64.(c)
    sv, U = coneig(sqrt.(c), sqrt.(conj.(c)), a, conj.(a))
    Um = U[:, 1:M]
    return balanced_truncation_sub(a, c, Um)
end
