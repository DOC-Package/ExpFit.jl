abstract type AbstractExpFit <: Function end
struct ExponentialFitting <: AbstractExpFit 
    expon::Vector{ComplexF64}
    coeff::Vector{ComplexF64}
end
struct ESPRIT <: AbstractExpFit end
struct MPencil <: AbstractExpFit end
struct Prony <: AbstractExpFit end

function (ef::ExponentialFitting)(t::Real)
    sum(ef.coeff .* exp.(-ef.expon .* t))
end

function expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; alg::AbstractExpFit=ESPRIT()) :: ExponentialFitting
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, nsamples, eps)
    elseif alg isa MPencil
        return matrix_pencil(func, tmin, tmax, nsamples, eps)
    elseif alg isa Prony
        return prony(func, tmin, tmax, nsamples, eps)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

struct ExponentialReduction <: AbstractExpFit
    expon::Vector{ComplexF64}
    coeff::Vector{ComplexF64}
end
struct BalancedTruncation <: AbstractExpFit end

function (ef::ExponentialReduction)(t::Real)
    sum(ef.coeff .* exp.(-ef.expon .* t))
end


function expred(exponent_init::AbstractVector{T}, coeff_init::AbstractVector{T}, eps::Real) where T<:ComplexF64
   
end