# ChebyJev
ChebyJev is the 
[OSCAR](https://www.oscar-system.org/)-based 
Julia version of the Maple package 
[GeneralizedChebyshev](https://github.com/TobiasMetzlaff/GeneralizedChebyshev). 

For Windows users, 
[WSL](https://learn.microsoft.com/en-us/windows/wsl/install) 
is required in order to access the Oscar functionalities. 

To install ChebyJev in Julia, run
```Julia
using Pkg
Pkg.add(url="https://github.com/TobiasMetzlaff/ChebyJev")
```
This will include the SDP instances from the example folder.

If you just want to use the functions package, download the ChebyJev.jl and run
```Julia
include("FILEPATH")
```
In WSL, use "/mnt/FILEPATH". 
Remember to add "ChebyJev." before any function call. 

## Dependencies
- [OSCAR](https://www.oscar-system.org/) (for algebraic operations with Chebyshev polynomials)
- [JuMP](https://github.com/jump-dev/JuMP.jl) (for semidefinite programming)

## Functionalities

### 1. Basics
Define the Lie type of the irreducible root system and its rank (A_n, B_n, C_n, D_n, E_6/7/8, F_4, G_2)
```Julia
using ChebyJev;
Type,n = :C,2;
```
Display base, coroots, positive roots, fundamental weights, highest root and more over Q (in the standard Euclidean basis) as depicted in Bourbaki, Groupes et Algebres de Lie, Ch. IV - VI, Planches I - IX
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
given by all z in R^n with H(z) a psd Hermite matrix polynomial with Hankel structure.
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
the chromatic number is bounded from below by the spectral bound 1-1/min f where f is a multiplicative invariant. 
Compute an order d semidefinite relaxation bound for the graph 
R^n with avoided set the Voronoi cell of the root lattice scaled by a factor r 
```Julia
using MosekTools
solver = Mosek
r=4
chromatic_euclidean_voronoi_sdp_bound(Type,n,d,r,solver)
```
or generate the SDP data in sdpa format to solve with another software
```Julia
folder = "FOLDERPATH"
chromatic_euclidean_voronoi_sdp_data(Type,n,d,r,folder)
```
Compute an order d semidefinite relaxation bound for the 
root lattice with avoided set the strict Voronoi vectors 
```Julia
using MosekTools
solver = Mosek
r=4
chromatic_lattice_sdp_bound(Type,n,d,r,solver)
```
or generate the SDP data in sdpa format to solve with another software
```Julia
folder = "FOLDERPATH"
chromatic_lattice_sdp_data(Type,n,d,r,folder)
```

## References
[1] [Optimization of trigonometric polynomials with crystallographic symmetry and spectral bounds for set avoiding graphs ](https://arxiv.org/abs/2303.09487)

[2] [Orbit Spaces of Weyl Groups Acting on Compact Tori: A Unified and Explicit Polynomial Description](https://arxiv.org/abs/2203.13152)

[3] [On symmetry adapted bases in trigonometric optimization](https://arxiv.org/abs/2310.05519)

[4] [Additive and Multiplicative Coinvariant Spaces of Weyl Groups in the Light of Harmonics and Graded Transfer](https://arxiv.org/abs/2412.17099)

## Contact
[Tobias Metzlaff](https://tobiasmetzlaff.com/): math"at"tobiasmetzlaff"dot"com
