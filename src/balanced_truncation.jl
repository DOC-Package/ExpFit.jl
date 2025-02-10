"""
    balanced_truncation_sub(a::Vector{ComplexF64}, c::Vector{ComplexF64}, Um::AbstractMatrix{ComplexF64})

Sub function for balanced truncation.
"""
function balanced_truncation_sub(a::Vector{T}, c::Vector{T}, Um::AbstractMatrix{ComplexF64}) where T<:Number

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
    balanced_truncation(a::Vector{ComplexF64}, c::Vector{ComplexF64}, eps::Float64)

Given the exponents `a` and coefficients `c` of an exponential fitting, compute a new exponential fitting with a reduced number of terms.
"""
function balanced_truncation(a::Vector{T}, c::Vector{T}, eps::Float64) where T<:Number
    a = ComplexF64.(a); c = ComplexF64.(c)
    sv, U = coneig(sqrt.(c), sqrt.(conj.(c)), a, conj.(a))
    M = findfirst(i -> 2 * sum(sv[i+1:end]) < eps, 1:length(sv))
    Um = U[:, 1:M]
    return balanced_truncation_sub(a, c, Um)
end

function balanced_truncation(a::Vector{T}, c::Vector{T}, M::Int) where T<:Number
    a = ComplexF64.(a); c = ComplexF64.(c)
    sv, U = coneig(sqrt.(c), sqrt.(conj.(c)), a, conj.(a))
    Um = U[:, 1:M]
    return balanced_truncation_sub(a, c, Um)
end
