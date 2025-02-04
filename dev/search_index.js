var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ExpFit","category":"page"},{"location":"#ExpFit","page":"Home","title":"ExpFit","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ExpFit.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ExpFit]","category":"page"},{"location":"#ExpFit.esprit-Tuple{AbstractVector{<:ComplexF64}, Real, Real}","page":"Home","title":"ExpFit.esprit","text":"esprit(hk, dt, eps; p=nothing) -> (exponent, coeff)\n\nPerform the ESPRIT algorithm using discrete data hk and the sampling interval dt.\n\nhk : Discrete data (vector of complex numbers).\ndt : Sampling interval.\neps: Threshold for singular value determination.\np  : Number of rows for the Hankel matrix (optional).\n\nReturns: A tuple containing the estimated exponents and coefficients. Note: The solution of the Vandermonde system is delegated to the function solve_vandermonde, which is assumed to be implemented elsewhere.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.esprit-Tuple{Function, Real, Real, Int64, Real}","page":"Home","title":"ExpFit.esprit","text":"esprit(func, tmin, tmax, N, eps; p=nothing) -> (exponent, coeff)\n\nSample the function func over the interval [tmin, tmax] to create discrete data, and then perform the ESPRIT algorithm.\n\nfunc : The function to be sampled (should return a real or complex number).\ntmin, tmax : Endpoints of the sampling interval.\nN    : Number of samples (this implementation generates 2N data points).\neps  : Threshold for singular value determination.\np    : Number of rows for the Hankel matrix (optional).\n\nReturns: A tuple containing the estimated exponents and coefficients.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.esprit_sub-Tuple{AbstractVector{<:ComplexF64}, Real}","page":"Home","title":"ExpFit.esprit_sub","text":"esprit_sub(hk, eps; p=nothing)\n\nEstimate the eigenvalues γ using the ESPRIT subspace method from the Hankel matrix constructed from hk.\n\nhk : A vector of complex data.\neps: Threshold used to distinguish noise. Singular values smaller than eps times the largest singular value are considered noise.\np  : Number of rows for the Hankel matrix (optional).\n\nReturns: A vector of eigenvalues γ.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.hankel_matrix-Tuple{AbstractVector{<:ComplexF64}}","page":"Home","title":"ExpFit.hankel_matrix","text":"hankel_matrix(hk; q=nothing)\n\nConstruct a Hankel matrix from the discrete data hk (a vector of complex numbers).\n\nhk : Data vector (of length K).\nq  : Number of columns for the Hankel matrix (optional).         If not provided, the default is:          - For K = 2N (even), q is set to N+1.          - For K = 2N-1 (odd),  q is set to N.        The allowed values are:          - For K = 2N: 2 ≤ q ≤ N+1.          - For K = 2N-1: 2 ≤ q ≤ N.\n\nReturns: A p × q Hankel matrix H such that H[i,j] = hk[i+j-1],          where p is computed as p = K + 1 - q.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.matrix_pencil-Tuple{AbstractVector{<:ComplexF64}, Real, Real}","page":"Home","title":"ExpFit.matrix_pencil","text":"matrix_pencil(hk, tmin, tmax, L, epsin) -> (exponent, coeff)\n\nPerform the Matrix Pencil method using discrete data hk and the sampling interval defined by [tmin, tmax].\n\nhk    : Discrete data (a vector of complex numbers), expected to be of length 2N.\ntmin  : The start of the sampling interval.\ntmax  : The end of the sampling interval.\nL     : Parameter for the Hankel matrix construction (the matrix will be of size (2N-L)×(L+1)).\nepsin : Threshold for determining the model order in the Matrix Pencil method.\n\nReturns: A tuple (exponent, coeff) where:\n\nexponent is a vector of estimated exponents (e.g. growth/decay rates).\ncoeff    is a vector of estimated coefficients.\n\nNote: The Vandermonde system is solved via least squares.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.matrix_pencil-Tuple{Function, Real, Real, Int64, Real}","page":"Home","title":"ExpFit.matrix_pencil","text":"matrix_pencil(func, tmin, tmax, K, L, epsin) -> (exponent, coeff)\n\nSample the function func over the interval [tmin, tmax] at K equally spaced points, and then perform the Matrix Pencil method.\n\nfunc  : The function to be sampled. It should accept a real argument t and return a ComplexF64.\ntmin  : The start of the sampling interval.\ntmax  : The end of the sampling interval.\nK     : Number of samples to generate (must be 2N, i.e. even).\nL     : Parameter for the Hankel matrix construction.\nepsin : Threshold for determining the model order.\n\nReturns: A tuple (exponent, coeff) containing the estimated exponents and coefficients.\n\n\n\n\n\n","category":"method"},{"location":"#ExpFit.matrix_pencil_sub-Tuple{AbstractVector{<:ComplexF64}, Int64, Real}","page":"Home","title":"ExpFit.matrix_pencil_sub","text":"matrix_pencil_sub(hk, L, epsin) -> (gamm, M)\n\nEstimate the eigenvalues γ using the Matrix Pencil method from the Hankel matrix constructed from the discrete data hk.\n\nhk    : A vector of complex samples. It is assumed that the length of hk is 2N.\nL     : An integer parameter for constructing the Hankel matrix. The Hankel matrix           is built with dimensions (2N - L) × (L + 1).\nepsin : A threshold used for model order determination. Specifically, the model order           M is determined as the smallest index m (with 1 ≤ m ≤ L) such that            |R[m+1, m+1]|/|R[1,1]| < epsin. In that case, M is set to m+1;           if no such index is found, M is set to 1.\n\nReturns: A tuple (gamm, M) where:\n\ngamm is a vector of eigenvalues obtained from the pencil.\nM is the detected model order.\n\n\n\n\n\n","category":"method"}]
}
