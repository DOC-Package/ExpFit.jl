struct ExponentialFitting <: Function 
    expon::Vector{ComplexF64}
    coeff::Vector{ComplexF64}
end

function expfit(func::Function, tmin::Real, tmax::Real, K::Int, eps::Real; method::ExponentialFitting)
    if method == ESPRIT
        return esprit(func, tmin, tmax, K, eps)
    elseif method == MatrixPencil
        return matrix_pencil(func, tmin, tmax, K, eps)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

function (ef::ExponentialFitting)(t::Real)
    sum(ef.coeff .* exp.(-ef.expon .* t))
end