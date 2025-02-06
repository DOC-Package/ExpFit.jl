struct ExponentialFitting <: Function 
    expon::Vector{ComplexF64}
    coeff::Vector{ComplexF64}
end

function (ef::ExponentialFitting)(t::Real)
    sum(ef.coeff .* exp.(-ef.expon .* t))
end

function expfit(func::Function, tmin::Real, tmax::Real, eps::Real; nsamples::Int=500, Smethod::ExponentialFitting=ESPRIT)
    if method == ESPRIT
        return esprit(func, tmin, tmax, K, eps)
    elseif method == MatrixPencil
        return matrix_pencil(func, tmin, tmax, K, eps)
    elseif method == Prony
        return prony(func, tmin, tmax, K, eps)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end




struct ExponentialReduction <: Function 
    expon::Vector{ComplexF64}
    coeff::Vector{ComplexF64}
end

function (ef::ExponentialReduction)(t::Real)
    sum(ef.coeff .* exp.(-ef.expon .* t))
end


function expred(exponent_init::AbstractVector{T}, coeff_init::AbstractVector{T}, eps::Real) where T<:ComplexF64
   
end