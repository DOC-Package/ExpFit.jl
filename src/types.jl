abstract type AbstractExpFit <: Function end

"""
    Exponentials <: AbstractExpFit

A type representing a sum of exponentials of the form

    f(t) = ∑ cᵢ exp(-aᵢ(t-t₀))

where `aᵢ` and `cᵢ` are the exponents and coefficients of the sum, respectively.

# Fields
- `expon::AbstractVector{<:Number}`: The exponents `aᵢ` of the sum of exponentials.
- `coeff::AbstractVector{<:Number}`: The coefficients `cᵢ` of the sum of exponentials.
"""
struct Exponentials <: AbstractExpFit 
    expon::AbstractVector{<:Number}
    coeff::AbstractVector{<:Number}
end

function (ef::Exponentials)(t::Real, t0::Real=0.0)
    sum(ef.coeff .* exp.(-ef.expon .* (t-t0)))
end

struct ESPRIT <: AbstractExpFit end
struct Pencil <: AbstractExpFit end
struct Prony <: AbstractExpFit end
struct ESPIRA1 <: AbstractExpFit end
struct ESPIRA2 <: AbstractExpFit end
struct FastESPRIT <: AbstractExpFit end
struct BalancedTruncation <: AbstractExpFit end


"""
    expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `func` using the algorithm `alg`.

# Arguments
- `func::Function`: The function to estimate the exponents and coefficients.
- `tmin::Real`: The initial time.
- `tmax::Real`: The final time.
- `nsamples::Int`: The number of samples to take in the interval [tmin, tmax].
- `eps::Real`: The tolerance for the algorithm.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, nsamples, eps)
    elseif alg isa Pencil
        return matrix_pencil(func, tmin, tmax, nsamples, eps)
    elseif alg isa Prony
        return prony(func, tmin, tmax, nsamples, eps)
    elseif alg isa ESPIRA1
        return espira1(func, tmin, tmax, nsamples, eps)
    elseif alg isa ESPIRA2
        return espira2(func, tmin, tmax, nsamples, eps)
    elseif alg isa FastESPRIT
        return fast_esprit(func, tmin, tmax, nsamples, eps)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

"""
    expfit(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `func` using the algorithm `alg`.

# Arguments
- `func::Function`: The function to estimate the exponents and coefficients.
- `tmin::Real`: The initial time.
- `tmax::Real`: The final time.
- `dt::Real`: The sampling interval.
- `eps::Real`: The tolerance for the algorithm.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(func::Function, tmin::Real, tmax::Real, dt::Real, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, dt, eps)
    elseif alg isa Pencil
        return matrix_pencil(func, tmin, tmax, dt, eps)
    elseif alg isa Prony
        return prony(func, tmin, tmax, dt, eps)
    elseif alg isa ESPIRA1
        return espira1(func, tmin, tmax, dt, eps)
    elseif alg isa ESPIRA2
        return espira2(func, tmin, tmax, dt, eps)
    elseif alg isa FastESPRIT
        return fast_esprit(func, tmin, tmax, dt, eps)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

"""
    expfit(f::AbstractVector{<:Number}, dt::Real, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `f` using the algorithm `alg`.

# Arguments
- `f::AbstractVector{<:Number}`: The samples of the function.
- `dt::Real`: The sampling interval.
- `eps::Real`: The tolerance for the algorithm.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(func::Function, tmin::Real, tmax::Real, nsamples::Int, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, nsamples, M)
    elseif alg isa Pencil
        return matrix_pencil(func, tmin, tmax, nsamples, M)
    elseif alg isa Prony
        return prony(func, tmin, tmax, nsamples, M)
    elseif alg isa ESPIRA1
        return espira1(func, tmin, tmax, nsamples, M)  
    elseif alg isa ESPIRA2
        return espira2(func, tmin, tmax, nsamples, M)
    elseif alg isa FastESPRIT
        return fast_esprit(func, tmin, tmax, nsamples, M)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

"""
    expfit(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `func` using the algorithm `alg`.

# Arguments
- `func::Function`: The function to estimate the exponents and coefficients.
- `tmin::Real`: The initial time.
- `tmax::Real`: The final time.
- `dt::Real`: The sampling interval.
- `M::Int`: The model order.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(func::Function, tmin::Real, tmax::Real, dt::Real, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(func, tmin, tmax, dt, M)
    elseif alg isa Pencil
        return matrix_pencil(func, tmin, tmax, dt, M)
    elseif alg isa Prony
        return prony(func, tmin, tmax, dt, M)
    elseif alg isa ESPIRA1
        return espira1(func, tmin, tmax, dt, M)  
    elseif alg isa ESPIRA2
        return espira2(func, tmin, tmax, dt, M)
    elseif alg isa FastESPRIT
        return fast_esprit(func, tmin, tmax, dt, M)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

"""
    expfit(f::AbstractVector{<:Number}, dt::Real, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `f` using the algorithm `alg`.

# Arguments
- `f::AbstractVector{<:Number}`: The samples of the function.
- `dt::Real`: The sampling interval.
- `eps::Real`: The tolerance for the algorithm.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(f::AbstractVector{<:Number}, dt::Real, eps::Real; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(f, dt, eps)
    elseif alg isa Pencil
        return matrix_pencil(f, dt, eps)
    elseif alg isa Prony
        return prony(f, dt, eps)
    elseif alg isa ESPIRA1
        return espira1(f, dt, eps)  
    elseif alg isa ESPIRA2
        return espira2(f, dt, eps)
    elseif alg isa FastESPRIT
        return fast_esprit(f, dt, eps)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end

"""
    expfit(f::AbstractVector{<:Number}, dt::Real, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials

Estimate the exponents and coefficients of the sum of exponentials for the function `f` using the algorithm `alg`.

# Arguments
- `f::AbstractVector{<:Number}`: The samples of the function.
- `dt::Real`: The sampling interval.
- `M::Int`: The model order.
- `alg::AbstractExpFit=ESPRIT()`: The algorithm to use for the estimation.
"""
function expfit(f::AbstractVector{<:Number}, M::Int; alg::AbstractExpFit=ESPRIT()) :: Exponentials
    if alg isa ESPRIT
        return esprit(f, dt, M)
    elseif alg isa Pencil
        return matrix_pencil(f, dt, M)
    elseif alg isa Prony
        return prony(f, dt, M)
    elseif alg isa ESPIRA1
        return espira1(f, dt, M)  
    elseif alg isa ESPIRA2
        return espira2(f, dt, M)
    elseif alg isa FastESPRIT
        return fast_esprit(f, dt, M)
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end


"""
    expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, eps::Real) :: Exponentials

Finding a new set of exponents and coefficients for the sum of exponentials using the balanced truncation method.

# Arguments
- `a::AbstractVector{<:Number}`: The exponents `aᵢ` of the sum of exponentials.
- `c::AbstractVector{<:Number}`: The coefficients `cᵢ` of the sum of exponentials.
- `eps::Real`: The tolerance for the balanced truncation method.
"""
function expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, eps::Real)
    return balanced_truncation(a, c, eps)
end

"""
    expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, M::Int) :: Exponentials

Finding a new set of exponents and coefficients for the sum of exponentials using the balanced truncation method.

# Arguments
- `a::AbstractVector{<:Number}`: The exponents `aᵢ` of the sum of exponentials.
- `c::AbstractVector{<:Number}`: The coefficients `cᵢ` of the sum of exponentials.
- `M::Int`: The model order.
"""
function expred(a::AbstractVector{<:Number}, c::AbstractVector{<:Number}, M::Int)
    return balanced_truncation(a, c, M)
end
