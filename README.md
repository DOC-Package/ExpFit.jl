# ExpFit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://doc-package.github.io/ExpFit.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://doc-package.github.io/ExpFit.jl/dev/)
[![Build Status](https://github.com/DOC-Package/ExpFit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DOC-Package/ExpFit.jl/actions/workflows/CI.yml?query=branch%3Amain)

Exponential Fitting Package

This package contains algorithms for fitting functions or discrete data with a sum of exponentials and for reducing the order of a sum of exponentials.

## Installation

```julia
    julia> using Pkg

    julia> Pkg.add(; url="https://github.com/DOC-package/ExpFit.jl")
```

## Simple usage

As an example, we consider approximating a Bessel function.  It is necessary to specify the range of approximation `[0.0,10.0]`, the number of sample points `N=100`, and the tolerance `tol=1e-2`.  Here is the script.

```julia
   julia> using ExpFit

   julia> using SpecialFunctions

   julia> f = t -> besselj(0,t) + 1.0im*besselj(1,t)
   #1 (generic function with 1 method)
   
   julia> ef = expfit(f, 0.0, 10.0, 100, 1e-2)
   (::Exponentials) (generic function with 2 methods)
```

The obtained exponents and coefficients are contained in `ef`
```julia
   julia> println("Exponents: ", ef.expon)
   Exponents: ComplexF64[0.6452980998429971 - 0.45960031394319034im, 0.07441134082707875 - 0.9779354644348793im, 0.4140265306909683 + 0.7881686630201671im]
   
   julia> println("Coefficients: ", ef.coeff)
   Coefficients: ComplexF64[0.5606982418231227 + 0.11793524228553061im, 0.43232179005370963 - 0.2865070855029491im, 0.013151533160065482 + 0.16322208134704955im]
```

As another use case, you can directly input equally spaced discrete sample data.  In this case, a time interval is required.
```julia
    julia> t = range(0.0, 10.0, length=100)
    0.0:0.10101010101010101:10.0

    julia> fv = f.(t)
    100-element Vector{ComplexF64}:
                      1.0 + 0.0im
       0.9974508660068557 + 0.05044066474846497im
       0.9898229555172593 + 0.10049567146911764im
                          ⋮
     -0.24027191897182099 + 0.06861955505415944im
      -0.2459357644513483 + 0.04347274616886144im

    julia> ef = expfit(fv, t[2]-t[1], 1e-2)
    (::Exponentials) (generic function with 2 methods)
```

## Documentation

For more details, please refer to [the documentation](https://doc-package.github.io/ExpFit.jl/dev/).

## Citation

The algorithms in this package are summerized in
- H. Takahashi, S. Rudge, C. Kaspar, M. Thoss, and R. Borrelli, J. Chem. Phys. 160, 204105 (2024). (https://doi.org/10.1063/5.0209348) 

For citations to the original papers on the algorithms, please check the references in the paper. 

If you find the package useful in your research, we would be grateful if you could cite our publication.
Here are the bibtex entries:
```bib
@article{TakahashiEtAl2024TJoCP,
  title = {High Accuracy Exponential Decomposition of Bath Correlation Functions for Arbitrary and Structured Spectral Densities: {{Emerging}} Methodologies and New Approaches},
  shorttitle = {High Accuracy Exponential Decomposition of Bath Correlation Functions for Arbitrary and Structured Spectral Densities},
  author = {Takahashi, Hideaki and Rudge, Samuel and Kaspar, Christoph and Thoss, Michael and Borrelli, Raffaele},
  year = {2024},
  month = may,
  journal = {The Journal of Chemical Physics},
  volume = {160},
  number = {20},
  pages = {204105},
  issn = {0021-9606},
  doi = {10.1063/5.0209348},
}

```




