using Random

function generate_exponent_coefficient_pairs(n::Int)
    Random.seed!(1234)
    exponents = Vector{ComplexF64}(undef, n)
    coefficients = Vector{ComplexF64}(undef, n)
    for i in 1:n
        re_exp = rand() * 9.9 + 0.1
        im_exp = randn()  
        exponents[i] = re_exp + im_exp*im

        re_coeff = randn()
        im_coeff = randn()
        coefficients[i] = re_coeff + im_coeff*im
    end
    return exponents, coefficients
end
