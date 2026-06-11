# ChebyJev

ChebyJev is an 
[OSCAR](https://www.oscar-system.org/)-based 
Julia version of the Maple package 
[GeneralizedChebyshev](https://github.com/TobiasMetzlaff/GeneralizedChebyshev). 

While generalized Chebyshev polynomials have wonderful properties for symbolical and numerical computation, 
their definition is based on the theory of root systems. 
This package aims to enable efficient computations for this context, 
in particular, for Weyl groups actions, multiplicative (co-)invariants and symmetry reduction. 



## Installation

For Windows users, 
[WSL](https://learn.microsoft.com/en-us/windows/wsl/install) 
is required in order to access the Oscar functionalities until Oscar becomes available for Windows. 

To install ChebyJev in Julia, run
```Julia
using Pkg
Pkg.add(url="https://github.com/TobiasMetzlaff/ChebyJev")
```
This will include the SDP instances from the example folder.
If you only want to use the functions, 
download the file ChebyJev.jl, run
```Julia
include("FILEPATH")
```
and add "ChebyJev." before any function call. 
In WSL, use "/mnt/FILEPATH". 



## Dependencies
- [OSCAR](https://www.oscar-system.org/) (for algebraic operations with Chebyshev polynomials)
- [JuMP](https://github.com/jump-dev/JuMP.jl) (for semidefinite programming)



## Credit

Algorithmic theory for Chebyshev polynomials and Weyl group orbits has been developed by 
Evelyne Hubert and Michael Singer in [6] and [7]. 
I am particularly grateful to Evelyne Hubert for suggestions, corrections and explanations concerning the implementation. 



## References
[1] [Optimization of trigonometric polynomials with crystallographic symmetry and spectral bounds for set avoiding graphs ](https://arxiv.org/abs/2303.09487)

[2] [Orbit Spaces of Weyl Groups Acting on Compact Tori: A Unified and Explicit Polynomial Description](https://arxiv.org/abs/2203.13152)

[3] [On symmetry adapted bases in trigonometric optimization](https://arxiv.org/abs/2310.05519)

[4] [Additive and Multiplicative Coinvariant Spaces of Weyl Groups in the Light of Harmonics and Graded Transfer](https://arxiv.org/abs/2412.17099)

[5] [Groupes et algebres de Lie, Ch. IV - VI](https://link.springer.com/book/10.1007/978-3-540-34491-9)

[6] [Sparse Interpolation in Terms of Multivariate Chebyshev Polynomials](https://link.springer.com/article/10.1007/s10208-021-09535-7)

[7] [Weyl group actions on minuscule weights in Maple](https://singer.math.ncsu.edu/papers/minuscule/)



## Contact
[Tobias Metzlaff](https://tobiasmetzlaff.com/): math"at"tobiasmetzlaff"dot"com / tobias"dot"metzlaff"at"bimsa"dot"cn



## Functionalities

### 1. Basics
Define the Lie type of the irreducible root system and its rank (A_n, B_n, C_n, D_n, E_6/7/8, F_4, G_2)
```Julia
using ChebyJev;
Type,n = :C,2;
```
Display base, coroots, positive roots, fundamental weights, highest root and more over Q (in the standard Euclidean basis) as depicted in Plates I - IX of [5]
```Julia
B = qbasematrix(Type,n)
typeof(B)
```
or Z (in the basis of fundamental weights)
```Julia
B = zbasematrix(Type,n)
typeof(B)
```

### 2. Weyl Group Action
Display the order of the Weyl group and its generating reflections over Q or Z
```Julia
weylgrouporder(Type,n)
qweylgroupgen(Type,n)
zweylgroupgen(Type,n)
```
Compute the orbit of a weight under the Weyl group
```Julia
M = qweightmatrix(Type,n)
weight = M[:,1]
orb1 = orbit(Type,n,weight)
alpha = [1,0]
orb2 = orbit(Type,alpha)
M*orb2   # should be orb1
```
Compute the unique representative of an orbit in the fundamental Weyl chamber
```Julia
M = qweightmatrix(Type,n)
weight = -2*M[:,1] + 1*M[:,2]
rep1 = pull(Type,n,weight)
alpha = [-2,1]
rep2 = pull(Type,alpha)
M*rep2   # should be rep1
```

### 3. Chebyshev Polynomials
Compute the generalized Chebyshev polynomial of the first kind associated to a weight
```Julia
M = qweightmatrix(Type,n)
weight = -2*M[:,1] + 1*M[:,2]
F = chebyshevpolynomial(Type,n,weight)
alpha = [-2,1]
chebyshevpolynomial(Type,alpha)  # should be F
parent(F)   # Multivariate polynomial ring in n variables z1, ..., zn over rational field
```

### 4. Multiplicative Invariants
Display the fundamental invariants of the muliplicative Weyl group action on Laurent polynomials
```Julia
theta = fundamentalinvariant(Type,n)
f = theta[1]^2*theta[2] + theta[1]
parent(f)   # Multivariate polynomial ring in n variables x1, .., xn, y1, ..., yn over rational field
```
Test if a Laurent polynomial is a multiplicative invariant
```Julia
is_multiplicative_invariant(Type,f)
```
Rewrite a multiplicative invariant in the fundamental invariant polynomial basis
```Julia
invariantrewrite_z(Type,f)
```
or basis of Chebyshev polynomials
```Julia
invariantrewrite_T(Type,f)
```

### 5. T-Orbit Spaces
The orbit space of the multiplicative Weyl group action is a compact basic semi algebraic set, 
given by all z in R^n with H(z) a psd Hermite matrix polynomial with Hankel structure [2].
Display the entry H_ij of H in the fundamental invariant polynomial basis
```Julia
i,j = 1,2
hermite_entry_z(Type,n,i+j)
```
or basis of Chebyshev polynomials
```Julia
hermite_entry_T(Type,n,i+j)
```
Display the whole matrix
```Julia
hermite_matrix_z(Type,n)
hermite_matrix_T(Type,n)
```

### 6. Polynomial Optimization
Compute the moment matrix in the basis of Chebyshev polynomials up to order d
```Julia
d=n+1
tmomentmatrix(Type,n,d)
```
Compute the localized moment matrix for optimization over the orbit space
```Julia
d=n+1
tlocalizedmomentmatrix(Type,n,d)
```
For certain set avoiding graphs with Weyl group symmetry, 
the chromatic number is bounded from below by the spectral bound 1-1/min f where f is a multiplicative invariant [1]. 
Compute an order d semidefinite relaxation bound for 
R^n with avoided set the Voronoi cell of the root lattice scaled by a factor r 
```Julia
using MosekTools
solver = Mosek
r=4
chromatic_euclidean_voronoi_sdp_bound(Type,n,d,r,solver)
```
or generate the SDP data in sdpa format directly 
```Julia
folder = "FOLDERPATH"
chromatic_euclidean_voronoi_sdp_data(Type,n,d,r,folder)
```
In this case, obtain an approximate spectral bound with 1+1/optimum. 
The change of sign comes from how sdpa encodes SDPs. 
Remember to include "/mnt/" for the folder when using WSL. 

Compute an order d semidefinite relaxation bound for the 
root lattice with avoided set the strict Voronoi vectors 
```Julia
using MosekTools
solver = Mosek
r=4
chromatic_lattice_sdp_bound(Type,n,d,r,solver)
```
or generate the SDP data in sdpa format directly 
```Julia
folder = "FOLDERPATH"
chromatic_lattice_sdp_data(Type,n,d,r,folder)
```
