abstract type AbstractExpFit <: Function end
struct Exponentials <: AbstractExpFit 
    expon::AbstractVector{ComplexF64}
    coeff::AbstractVector{ComplexF64}
end
struct ESPRIT <: AbstractExpFit end
struct MPencil <: AbstractExpFit end
struct Prony <: AbstractExpFit end
struct BalancedTruncation <: AbstractExpFit end

function (ef::Exponentials)(t::Real, t0::Real=0.0)
    sum(ef.coeff .* exp.(-ef.expon .* (t-t0)))
end

function expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials
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

function expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, nsamples, M)
    elseif alg isa MPencil
        return matrix_pencil(func, tmin, tmax, nsamples, M)
    elseif alg isa Prony
        return prony(func, tmin, tmax, nsamples, M)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

function expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, eps::Real)
    return balanced_truncation(a, c, eps)
end

function expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, M::Int)
    return balanced_truncation(a, c, M)
end
