module ChebyJev

using Oscar

# ============================================================
# ============================================================
# Internal Functions
# ============================================================
# ============================================================

"""
    leftrightinverse(M)

Compute a left or right inverse depending on the shape of the matrix.

- If rows <= columns, returns a left inverse: (M*M')^-1 * M
- If rows >= columns, returns a right inverse: M'*(M*M')^-1
"""
function leftrightinverse(M)
    if nrows(M) <= ncols(M)
        G = M * transpose(M)
        return transpose(M) * inv(G)
    else
        G = transpose(M) * M
        return inv(G) * transpose(M)
    end
end

#############

"""
    coroot(r)

Return the coroot corresponding to the root `r`.

h = 2r/(r,r), where (r,r) is the Euclidean scalar product.
"""
function coroot(r::AbstractVector{<:QQFieldElem})
    scalar = sum(x^2 for x in r)
    return (QQ(2)//scalar) .* r
end

#############

"""
    esp(L::Vector, r::Int)

Compute the r-th elementary symmetric polynomial in the entries of L.

- L: a vector of elements (numbers or symbols)
- r: the order of the elementary symmetric polynomial
"""
function esp(L::Vector, r::Int)
    n = length(L)
    if r == 0
        return one(parent(L[1]))
    elseif r < 0 || r > n
        return zero(parent(L[1]))
    end
    e = [zero(parent(L[1])) for i in 1:r+1]
    e[1] = one(parent(L[1]))
    for a in L
        for k in min(r,n):-1:1
            e[k+1] += a*e[k]
        end
    end
    return e[r+1]
end

#############

# ============================================================
# ============================================================
# Exported Functions
# ============================================================
# ============================================================

export qbasematrix,
       zbasematrix,
       qpositiveroots,
       zpositiveroots,
       qweightmatrix,
       zweightmatrix,
       qhighestroot,
       zhighestroot,
       qweylgroupgen,
       zweylgroupgen,
       orbitcardinality,
       weyllength,
       steinbergweight,
       pull,
       orbit,
       chebyshev,
       fundamental_invariant,
       invariant_rewrite,
       moment_matrix_T,
       localized_moment_matrix_T,
       chromatic_sdp_bound,
       chromatic_sdp_data

# ============================================================
# 1. Basics
# ============================================================

function qtozcoords(Type::Symbol, n::Int, M)

    W = qweightmatrix(Type,n)
    L = leftrightinverse(W)

    A = L * M

    Z = Matrix{Int64}(undef, nrows(A), ncols(A))

    for i in 1:nrows(A)
        for j in 1:ncols(A)
            @assert denominator(A[i,j]) == 1
            Z[i,j] = Int64(A[i,j])
        end
    end

    return Z

end

#############

"""
    qbasematrix(Type::Symbol, n::Int)

Return the matrix whose columns are the simple roots of the
root system of type `Type` and rank `n`.

Each column is a simple root vector.
"""
function qbasematrix(Type::Symbol, n::Int)
    roots = []

    if Type == :A
        for i in 1:n
            push!(roots, [QQ(j == i ? 1 : (j == i+1 ? -1 : 0)) for j in 1:n+1])
        end

    elseif Type == :B
        for i in 1:n
            push!(roots, [QQ(j == i ? 1 : (j == i+1 ? -1 : 0)) for j in 1:n])
        end

    elseif Type == :C
        for i in 1:n-1
            push!(roots, [QQ(j == i ? 1 : (j == i+1 ? -1 : 0)) for j in 1:n])
        end
        push!(roots, [QQ(j == n ? 2 : 0) for j in 1:n])

    elseif Type == :D
        for i in 1:n-1
            push!(roots, [QQ(j == i ? 1 : (j == i+1 ? -1 : 0)) for j in 1:n])
        end
        push!(roots, [QQ(j == n ? 1 : (j == n-1 ? 1 : 0)) for j in 1:n])

    elseif Type == :E && 6 <= n <= 8

        push!(roots, [
            QQ(1)//2, QQ(-1)//2, QQ(-1)//2, QQ(-1)//2,
            QQ(-1)//2, QQ(-1)//2, QQ(-1)//2, QQ(1)//2
        ])

        push!(roots, [
            QQ(1), QQ(1), QQ(0), QQ(0),
            QQ(0), QQ(0), QQ(0), QQ(0)
        ])

        for i in 3:n
            push!(roots, [
                QQ(j == i - 1 ? 1 : (j == i - 2 ? -1 : 0))
                for j in 1:8
            ])
        end

    elseif Type == :F && n == 4
        push!(roots, [QQ(0),QQ(1),QQ(-1),QQ(0)])
        push!(roots, [QQ(0),QQ(0),QQ(1),QQ(-1)])
        push!(roots, [QQ(0),QQ(0),QQ(0),QQ(1)])
        push!(roots, [QQ(1)//2,QQ(-1)//2,QQ(-1)//2,QQ(-1)//2])

    elseif Type == :G && n == 2
        push!(roots, [QQ(1),QQ(-1),QQ(0)])
        push!(roots, [QQ(-2),QQ(1),QQ(1)])

    else
        error("Root system must be of simple Lie type")
    end

    return matrix(QQ, hcat(roots...))
end

"""
    zbasematrix(Type::Symbol, n::Int)

Return the simple roots in coordinates with respect to the basis of
fundamental weights.

The output is an Int64 matrix.
"""
function zbasematrix(Type::Symbol, n::Int)

    return qtozcoords(Type,n,qbasematrix(Type,n))

end

#############

"""
    qpositiveroots(Type::Symbol, n::Int)

Return a matrix whose columns are the positive roots of the
root system of type `Type` and rank `n`.
"""
function qpositiveroots(Type::Symbol, n::Int)
    B = qbasematrix(Type,n)
    roots = Vector{Vector{QQFieldElem}}()

    if Type == :A
        for i in 1:n
            for j in i+1:size(B,1)
                push!(roots, [QQ(k==i ? 1 : (k==j ? -1 : 0)) for k in 1:size(B,1)])
            end
        end

    elseif Type == :B
        for i in 1:n
            push!(roots, [QQ(k==i ? 1 : 0) for k in 1:n])
        end
        for i in 1:n-1
            for j in i+1:n
                push!(roots, [QQ(k==i ? 1 : (k==j ? -1 : 0)) for k in 1:n])
                push!(roots, [QQ(k==i ? 1 : (k==j ? 1 : 0)) for k in 1:n])
            end
        end

    elseif Type == :C
        for i in 1:n
            push!(roots, [QQ(k==i ? 2 : 0) for k in 1:n])
        end
        for i in 1:n-1
            for j in i+1:n
                push!(roots, [QQ(k==i ? 1 : (k==j ? -1 : 0)) for k in 1:n])
                push!(roots, [QQ(k==i ? 1 : (k==j ? 1 : 0)) for k in 1:n])
            end
        end

    elseif Type == :D
        for i in 1:n-1
            for j in i+1:n
                push!(roots, [QQ(k==i ? 1 : (k==j ? -1 : 0)) for k in 1:n])
                push!(roots, [QQ(k==i ? 1 : (k==j ? 1 : 0)) for k in 1:n])
            end
        end

    elseif Type == :F && n==4
        for i in 1:4
            push!(roots, [QQ(k==i ? 1 : 0) for k in 1:4])
        end
        for i in 1:3
            for j in i+1:4
                push!(roots, [QQ(k==i ? 1 : (k==j ? -1 : 0)) for k in 1:4])
                push!(roots, [QQ(k==i ? 1 : (k==j ? 1 : 0)) for k in 1:4])
            end
        end
        hs = [
            [1,1,1,1],[1,1,1,-1],[1,1,-1,1],[1,-1,1,1],
            [1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1],[1,-1,-1,-1]
        ]
        for v in hs
            push!(roots, [QQ(x)//2 for x in v])
        end

    elseif Type == :G && n==2
        push!(roots, collect(B[:,1]))
        push!(roots, collect(B[:,1]+B[:,2]))
        push!(roots, collect(2*B[:,1]+B[:,2]))
        push!(roots, collect(3*B[:,1]+B[:,2]))
        push!(roots, collect(3*B[:,1]+2*B[:,2]))
        push!(roots, collect(B[:,2]))
    else
        error("Root system must be of simple Lie type")
    end

    return matrix(QQ, hcat(roots...))
end

"""
    zpositiveroots(Type::Symbol, n::Int)

Return the positive roots in coordinates with respect to the basis of
fundamental weights.

The output is an Int64 matrix.
"""
function zpositiveroots(Type::Symbol, n::Int)

    return qtozcoords(Type,n,qpositiveroots(Type,n))

end

#############

"""
    qweightmatrix(Type::Symbol, n::Int)

Return the matrix whose columns are the fundamental weights of the
root system of type `Type` and rank `n`.

Columns correspond to fundamental weights. Returns a `Nemo.QQMatrix`
to be consistent with `qbasematrix`.
"""
function qweightmatrix(Type::Symbol, n::Int)

    B = qbasematrix(Type,n)

    # simple coroots as columns
    H = hcat([coroot(B[:,i]) for i in 1:n]...)

    # left inverse: columns = fundamental weights
    W_array = H * inv(transpose(H)*H)

    # convert to QQMatrix
    return matrix(QQ, W_array)
end

"""
    zweightmatrix(Type::Symbol, n::Int)

Return the fundamental weights in coordinates with respect to the basis
of fundamental weights.

The output is the Int64 identity matrix.
"""
function zweightmatrix(Type::Symbol, n::Int)

    return qtozcoords(Type,n,qweightmatrix(Type,n))

end

#############

"""
    qhighestroot(Type::Symbol, n::Int)

Return the highest root of the root system of type `Type` and rank `n`.
"""
function qhighestroot(Type::Symbol, n::Int)

    B = qbasematrix(Type,n)
    W = qweightmatrix(Type,n)

    if Type == :A

        return W[:,1] + W[:,n]

    elseif Type == :B

        return B[:,1] + 2 * sum(B[:,j] for j in 2:n)

    elseif Type == :C

        return 2 * W[:,1]

    elseif Type == :D

        return B[:,1] + 2 * sum(B[:,j] for j in 2:n-2) + B[:,n-1] + B[:,n]

    elseif Type == :E

        if n == 6
            return W[:,2]
        elseif n == 7
            return W[:,1]
        elseif n == 8
            return W[:,8]
        else
            error("Invalid rank for type E")
        end

    elseif Type == :F && n == 4

        return W[:,1]

    elseif Type == :G && n == 2

        return W[:,2]

    else

        error("Root system must be of simple Lie type")

    end

end

"""
    zhighestroot(Type::Symbol, n::Int)

Return the highest root in coordinates with respect to the basis of
fundamental weights.

The output is an Int64 vector.
"""
function zhighestroot(Type::Symbol, n::Int)

    v = qhighestroot(Type,n)
    M = matrix(QQ, length(v), 1, v)

    Z = qtozcoords(Type,n,M)

    return [Z[i,1] for i in 1:size(Z,1)]

end

#############

"""
    vertexfunddomcoefficient(Type::Symbol, n::Int)

Return the list of scalar divisors defining the vertices of the
fundamental domain.

The output consists of the coefficients of the highest root with
respect to the fundamental weights, followed by 1.
"""
function vertexfunddomcoefficient(Type::Symbol, n::Int)

    theta = qhighestroot(Type,n)
    W = qweightmatrix(Type,n)

    coeffs = [dot(theta, W[:,i]) for i in 1:ncols(W)]

    return vcat(coeffs, QQ(1))

end

#############

const CHEBYSHEV_LEVEL_CACHE =
    Dict{Tuple{Symbol,Int64,Int64},Vector{Vector{Int64}}}()

function chebyshevlevel(Type::Symbol, n::Int, l::Int)

    key = (Type,n,l)

    if haskey(CHEBYSHEV_LEVEL_CACHE,key)
        return CHEBYSHEV_LEVEL_CACHE[key]
    end

    Ffull = vertexfunddomcoefficient(Type,n)
    F = Ffull[1:n]

    levelscale = minimum(F)
    target = QQ(l) * levelscale

    out = Vector{Vector{Int64}}()
    v = zeros(Int64,n)

    function rec(i::Int, remaining)

        if i > n
            if remaining == 0
                push!(out, copy(v))
            end
            return
        end

        max_i = floor(Int, remaining / F[i])

        for a in 0:max_i
            v[i] = a
            rec(i+1, remaining - QQ(a)*F[i])
        end

        v[i] = 0

    end

    rec(1,target)

    CHEBYSHEV_LEVEL_CACHE[key] = out

    return out

end

function chebyshevbasis_upto(Type::Symbol, n::Int, d::Int)

    L = Vector{Vector{Int64}}()

    for l in 0:d
        append!(L, chebyshevlevel(Type,n,l))
    end

    return L

end

#############

"""
    invariantdegrees(Type::Symbol, n::Int)

Return a vector of the degrees of the fundamental invariants of
the Weyl group of the specified type and rank.
"""
function invariantdegrees(Type::Symbol, n::Int)

    if Type == :A
        return collect(2:n+1)

    elseif Type == :B || Type == :C
        return [2*i for i in 1:n]

    elseif Type == :D
        return vcat([2*i for i in 1:n-1], n)

    elseif Type == :E
        if n == 6
            return [2,5,6,8,9,12]
        elseif n == 7
            return [2,6,8,10,12,14,18]
        elseif n == 8
            return [2,8,12,14,18,20,24,30]
        else
            error("Invalid rank for type E")
        end

    elseif Type == :F && n == 4
        return [2,6,8,12]

    elseif Type == :G && n == 2
        return [2,6]

    else
        error("Root system must be of simple Lie type")
    end

end

#############

"""
    weylgrouporder(Type::Symbol, n::Int)

Return the order of the Weyl group of the specified type and rank.

Computed as the product of the invariant degrees.
"""
function weylgrouporder(Type::Symbol, n::Int)
    prod(invariantdegrees(Type,n))
end

#############

"""
    qweylgroupgen(Type::Symbol, n::Int)

Return a vector of generators of the Weyl group as orthogonal matrices
acting on the ambient space.

Each generator is the qreflection

    x -> x - dot(x,h)*a

where a is a simple root and h is its coroot.
"""
function qweylgroupgen(Type::Symbol, n::Int)

    B = qbasematrix(Type,n)

    gens = Vector{QQMatrix}()

    for i in 1:n

        a = B[:,i]
        h = coroot(a)

        m = length(a)

        A = matrix(QQ, m, 1, collect(a))
        H = matrix(QQ, 1, m, collect(h))

        S = identity_matrix(QQ,m) - A * H

        push!(gens, S)

    end

    return gens

end

#############

"""
    zweylgroupgen(Type::Symbol, n::Int)

Return a vector of generators of the Weyl group as integer matrices
acting on coordinates in the basis of fundamental weights.

The computation is done over QQ and converted to Int64 at the end.
"""
function zweylgroupgen(Type::Symbol, n::Int)

    W = qweightmatrix(Type,n)
    R = qweylgroupgen(Type,n)
    L = leftrightinverse(W)

    gens = Vector{Matrix{Int64}}()

    for s in R

        A = L * s * W

        M = Matrix{Int64}(undef, n, n)

        for i in 1:n
            for j in 1:n

                @assert denominator(A[i,j]) == 1

                M[i,j] = Int64(A[i,j])

            end
        end

        push!(gens, M)

    end

    return gens

end

#############



# ============================================================
# 2. Reflections and Weights
# ============================================================

"""
    weyllength(Type::Symbol, n::Int, s::QQMatrix)

Return the Weyl length of the Weyl group element `s`.

The input `s` is a matrix over QQ acting on the standard coordinates
used by the q-functions.
"""
function weyllength(Type::Symbol, n::Int, s::QQMatrix)

    P = qpositiveroots(Type,n)

    len = 0

    for i in 1:ncols(P)

        r = P[:,i]
        sr = s * r

        if any(sr == -P[:,j] for j in 1:ncols(P))
            len += 1
        end

    end

    return len

end

"""
    weyllength(Type::Symbol, B::Matrix{Int64})

Return the Weyl length of the Weyl group element `B`.

The input `B` is a matrix over Int64 acting on coordinates with respect
to the basis of fundamental weights.
"""
function weyllength(Type::Symbol, B::Matrix{Int64})

    n = size(B,1)

    P = zpositiveroots(Type,n)

    len = 0

    for i in 1:size(P,2)

        r = P[:,i]
        sr = B * r

        if any(sr == -P[:,j] for j in 1:size(P,2))
            len += 1
        end

    end

    return len

end

#############

"""
    qreflection(rho, omega)

Return the qreflection of the vector `omega` across the hyperplane
orthogonal to the root `rho`.
"""
function qreflection(rho::AbstractVector{<:QQFieldElem}, omega::AbstractVector{<:QQFieldElem})
    omega - dot(coroot(rho), omega) * rho
end

#############

"""
    orbit(Type::Symbol, n::Int, omega)

Return the orbit of `omega` under the Weyl group generated by
qreflections in the simple roots.

The output is a vector of vectors.
"""
function orbit(Type::Symbol, n::Int, omega::AbstractVector{<:QQFieldElem})

    B = qbasematrix(Type,n)

    orb = Set{Vector{QQFieldElem}}()
    stack = Vector{Vector{QQFieldElem}}()

    omega0 = collect(omega)

    push!(orb, omega0)
    push!(stack, omega0)

    while !isempty(stack)

        next_omega = pop!(stack)

        for i in 1:n

            new_omega = collect(qreflection(B[:,i], next_omega))

            if !(new_omega in orb)
                push!(orb, new_omega)
                push!(stack, new_omega)
            end

        end

    end

    return collect(orb)

end

function orbit(Type::Symbol, alpha::Vector{Int64})

    n = length(alpha)

    alphaQQ = [QQ(a) for a in alpha]

    W = qweightmatrix(Type,n)
    R = qweylgroupgen(Type,n)
    L = leftrightinverse(W)

    omega = W * matrix(QQ, n, 1, alphaQQ)

    orb = Set{Vector{QQFieldElem}}()
    stack = Vector{Vector{QQFieldElem}}()

    omega0 = [omega[i,1] for i in 1:nrows(omega)]

    push!(orb, omega0)
    push!(stack, omega0)

    while !isempty(stack)

        v = pop!(stack)
        vmat = matrix(QQ, length(v), 1, v)

        for s in R

            new_mat = s * vmat
            new_v = [new_mat[i,1] for i in 1:nrows(new_mat)]

            if !(new_v in orb)
                push!(orb, new_v)
                push!(stack, new_v)
            end

        end

    end

    out = Vector{Vector{Int64}}()

    for v in orb

        vmat = matrix(QQ, length(v), 1, v)
        amat = L * vmat

        for i in 1:nrows(amat)
            @assert denominator(amat[i,1]) == 1
        end

        push!(out, [Int64(amat[i,1]) for i in 1:nrows(amat)])

    end

    return out

end

#############

"""
    steinbergweight(Type::Symbol, n::Int, s::QQMatrix)

Compute the Steinberg weight of the Weyl group element `s`.

- `s` is a QQMatrix from qweylgroupgen.
- Returns a vector of Int64 in the fundamental weight basis.
"""
function steinbergweight(Type::Symbol, n::Int, s::QQMatrix)

    # fundamental weights as columns
    omega = [qweightmatrix(Type,n)[:,i] for i in 1:n]

    # generators of the Weyl group over QQ
    Gens = qweylgroupgen(Type,n)

    # select indices i where length decreases under s*Gens[i]
    L = [i for i in 1:n if weyllength(Type,n, s * Gens[i]) < weyllength(Type,n,s)]

    # sum of corresponding fundamental weights
    v = zeros(QQ, size(omega[1],1))
    for i in L
        v += omega[i]
    end

    return v

end

"""
    steinbergweight(Type::Symbol, B::Matrix{Int64})

Compute the Steinberg weight of the Weyl group element `B`.

The input `B` is a matrix over Int64 acting on coordinates with respect
to the basis of fundamental weights. The computation is done over QQ
and the output is converted back to Int64.
"""
function steinbergweight(Type::Symbol, B::Matrix{Int64})

    n = size(B,1)

    W = qweightmatrix(Type,n)
    M = leftrightinverse(W)
    BB = matrix(QQ, n, n, [QQ(B[i,j]) for i in 1:n for j in 1:n])
    s = W * BB * M
    v = steinbergweight(Type,n,s)
    v = M*v
    v = Vector{Int64}(round.(Int, v))

    return v

end

#############

"""
    pull(Type::Symbol, x::Vector{QQFieldElem}, n::Int)

Return the representative of the Weyl orbit of the ambient vector `x`
lying in the fundamental Weyl chamber.

The input and output are ambient vectors over QQ.
"""
function pull(Type::Symbol, n::Int, x::Vector{QQFieldElem})

    W = qweightmatrix(Type,n)
    L = leftrightinverse(W)
    B = qbasematrix(Type,n)

    mu = copy(x)

    while true

        mumat = matrix(QQ, length(mu), 1, mu)
        coords = L * mumat

        signs = [coords[i,1] < 0 ? -1 : (coords[i,1] > 0 ? 1 : 0) for i in 1:n]

        if !any(s == -1 for s in signs)
            return mu
        end

        for i in 1:n
            if signs[i] == -1
                mu = qreflection(B[:,i], mu)
            end
        end

    end

end

"""
    pull(Type::Symbol, alpha::Vector{Int64})

Return the representative of the Weyl orbit of `alpha` lying in the
fundamental Weyl chamber.

The input and output are coordinates with respect to the basis of
fundamental weights. The computation is done over QQ and the output is
converted to Int64.
"""
function pull(Type::Symbol, alpha::Vector{Int64})

    n = length(alpha)

    W = qweightmatrix(Type,n)
    L = leftrightinverse(W)

    alphaQQ = [QQ(a) for a in alpha]
    alphamat = matrix(QQ, n, 1, alphaQQ)

    mu = W * alphamat
    muvec = [mu[i,1] for i in 1:nrows(mu)]

    while true

        mumat = matrix(QQ, length(muvec), 1, muvec)
        coords = L * mumat

        signs = [coords[i,1] < 0 ? -1 : (coords[i,1] > 0 ? 1 : 0) for i in 1:n]

        if !any(s == -1 for s in signs)
            break
        end

        B = qbasematrix(Type,n)

        for i in 1:n
            if signs[i] == -1
                muvec = qreflection(B[:,i], muvec)
            end
        end

    end

    pulledmat = matrix(QQ, length(muvec), 1, muvec)
    coords = L * pulledmat

    out = Vector{Int64}()

    for i in 1:nrows(coords)
        @assert denominator(coords[i,1]) == 1
        push!(out, Int64(coords[i,1]))
    end

    return out

end

#############

function dynkin_edges(Type::Symbol, n::Int)

    B = qbasematrix(Type,n)
    H = hcat([coroot(B[:,i]) for i in 1:n]...)

    edges = Tuple{Int64,Int64,Int64}[]

    for i in 1:n
        for j in i+1:n

            cij = sum(H[k,i] * B[k,j] for k in 1:nrows(B))
            cji = sum(H[k,j] * B[k,i] for k in 1:nrows(B))

            p = cij * cji

            if p == 0
                continue
            elseif p == 1
                push!(edges, (i,j,3))
            elseif p == 2
                push!(edges, (i,j,4))
            elseif p == 3
                push!(edges, (i,j,6))
            else
                error("Invalid Dynkin edge")
            end

        end
    end

    return edges

end

function component_order(comp::Vector{Int64}, edges)

    r = length(comp)

    if r == 0
        return 1
    elseif r == 1
        return 2
    end

    cset = Set(comp)
    subedges = [(i,j,m) for (i,j,m) in edges if i in cset && j in cset]

    labels = [m for (i,j,m) in subedges]

    deg = Dict(i => 0 for i in comp)

    for (i,j,m) in subedges
        deg[i] += 1
        deg[j] += 1
    end

    if any(m -> m == 6, labels)

        if r == 2
            return weylgrouporder(:G,2)
        else
            error("Invalid component with edge label 6")
        end

    elseif any(m -> m == 4, labels)

        if r == 2
            return 8
        end

        four_edges = [(i,j) for (i,j,m) in subedges if m == 4]

        if length(four_edges) != 1
            error("Invalid component with edge label 4")
        end

        a, b = four_edges[1]

        if r == 4 && deg[a] == 2 && deg[b] == 2
            return weylgrouporder(:F,4)
        else
            return 2^r * factorial(r)
        end

    else

        maxdeg = maximum(collect(values(deg)))

        if maxdeg <= 2
            return factorial(r + 1)
        end

        branch = [i for i in comp if deg[i] == 3]

        if length(branch) != 1
            error("Invalid simply laced component")
        end

        b = branch[1]
        arms = Int64[]

        neighbors = Int64[]

        for (i,j,m) in subedges
            if i == b
                push!(neighbors, j)
            elseif j == b
                push!(neighbors, i)
            end
        end

        for u in neighbors

            prev = b
            cur = u
            len = 1

            while deg[cur] == 2

                next_vertices = Int64[]

                for (i,j,m) in subedges
                    if i == cur && j != prev
                        push!(next_vertices, j)
                    elseif j == cur && i != prev
                        push!(next_vertices, i)
                    end
                end

                prev = cur
                cur = next_vertices[1]
                len += 1

            end

            push!(arms, len)

        end

        sort!(arms)

        if arms[1] == 1 && arms[2] == 1
            return weylgrouporder(:D,r)
        elseif arms == [1,2,2]
            return weylgrouporder(:E,6)
        elseif arms == [1,2,3]
            return weylgrouporder(:E,7)
        elseif arms == [1,2,4]
            return weylgrouporder(:E,8)
        else
            error("Unknown simply laced component")
        end

    end

end

function parabolic_order(Type::Symbol, n::Int, alpha::Vector{Int64})

    zero_vertices = [i for i in 1:n if alpha[i] == 0]

    if isempty(zero_vertices)
        return 1
    end

    edges = dynkin_edges(Type,n)

    vset = Set(zero_vertices)
    adj = Dict(i => Int64[] for i in zero_vertices)

    for (i,j,m) in edges
        if i in vset && j in vset
            push!(adj[i], j)
            push!(adj[j], i)
        end
    end

    seen = Set{Int64}()
    orders = Int64[]

    for v in zero_vertices

        if v in seen
            continue
        end

        stack = [v]
        comp = Int64[]

        while !isempty(stack)

            x = pop!(stack)

            if x in seen
                continue
            end

            push!(seen, x)
            push!(comp, x)

            for y in adj[x]
                if !(y in seen)
                    push!(stack, y)
                end
            end

        end

        push!(orders, component_order(comp, edges))

    end

    return prod(orders)

end

function orbitcardinality(Type::Symbol, n::Int, omega::Vector{QQFieldElem})

    if all(x -> x == 0, omega)
        return 1
    end

    mu = pull(Type,n,omega)

    W = qweightmatrix(Type,n)
    L = leftrightinverse(W)

    mumat = matrix(QQ, length(mu), 1, mu)
    coords = L * mumat

    alpha = Int64[]

    for i in 1:nrows(coords)

        if coords[i,1] > 0
            push!(alpha, 1)
        elseif coords[i,1] == 0
            push!(alpha, 0)
        else
            error("pull did not return a dominant vector")
        end

    end

    Worder = weylgrouporder(Type,n)
    stabilizer_order = parabolic_order(Type,n,alpha)

    return div(Worder, stabilizer_order)

end

function orbitcardinality(Type::Symbol, beta::Vector{Int64})

    if all(x -> x == 0, beta)
        return 1
    end

    n = length(beta)

    W = qweightmatrix(Type,n)

    betamat = matrix(QQ, n, 1, [QQ(x) for x in beta])
    omega = W * betamat
    omegavec = [omega[i,1] for i in 1:nrows(omega)]

    return orbitcardinality(Type,n,omegavec)

end

#############

# ============================================================
# 3. Chebyshev Polynomials
# ============================================================

const CHEBYSHEV_CONTEXT = Dict{
    Tuple{Symbol,Int64}, 
    Tuple{QQMPolyRing, Vector{QQMPolyRingElem}, Dict{Tuple{Vararg{Int64}},QQMPolyRingElem}}
}()

const CHEB_MULT_CACHE = Dict{
    Tuple{Symbol,Tuple{Vararg{Int64}},Tuple{Vararg{Int64}}},
    Dict{Tuple{Vararg{Int64}},QQFieldElem}
}()

function chebyshevcontext(Type::Symbol, n::Int)

    key = (Type,n)

    if haskey(CHEBYSHEV_CONTEXT, key)
        return CHEBYSHEV_CONTEXT[key]
    end

    R, z = polynomial_ring(QQ, ["z$i" for i in 1:n])
    cache = Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem}()

    CHEBYSHEV_CONTEXT[key] = (R, collect(z), cache)

    return CHEBYSHEV_CONTEXT[key]

end

function cheb_mult(Type::Symbol,
                   alpha::Vector{Int64},
                   beta::Vector{Int64})

    key = (Type, Tuple(alpha), Tuple(beta))

    if haskey(CHEB_MULT_CACHE,key)
        return CHEB_MULT_CACHE[key]
    end

    orb = orbit(Type,beta)
    m = orbitcardinality(Type,beta)

    out = Dict{Tuple{Vararg{Int64}},QQFieldElem}()

    for gamma in orb

        delta = Tuple(pull(Type, alpha .+ gamma))

        out[delta] = get(out, delta, QQ(0)) + QQ(1)//QQ(m)

    end

    CHEB_MULT_CACHE[key] = out

    return out

end

#########

const TBASIS_ELEM = Dict{
    Tuple{Vararg{Int64}}, 
    QQFieldElem
}

function tbasis_key(Type::Symbol, alpha::Vector{Int64})
    return Tuple(pull(Type, alpha))
end

function tbasis_one(n::Int)
    return TBASIS_ELEM(Tuple(zeros(Int64,n)) => QQ(1))
end

function tbasis_term(Type::Symbol, alpha::Vector{Int64}, c = QQ(1))

    key = tbasis_key(Type, alpha)

    if c == 0
        return TBASIS_ELEM()
    end

    return TBASIS_ELEM(key => QQ(c))

end

function tbasis_add(a::TBASIS_ELEM, b::TBASIS_ELEM)

    out = TBASIS_ELEM()

    for (k,c) in a
        out[k] = get(out, k, QQ(0)) + c
    end

    for (k,c) in b
        out[k] = get(out, k, QQ(0)) + c
        if out[k] == 0
            delete!(out, k)
        end
    end

    return out

end

function tbasis_scale(c, a::TBASIS_ELEM)

    cq = QQ(c)

    if cq == 0
        return TBASIS_ELEM()
    end

    out = TBASIS_ELEM()

    for (k,d) in a
        cd = cq*d
        if cd != 0
            out[k] = cd
        end
    end

    return out

end

function tbasis_mul_basis(Type::Symbol, alpha::Tuple, beta::Tuple)

    return deepcopy(
        cheb_mult(Type, collect(alpha), collect(beta))
    )

end

function tbasis_mul(Type::Symbol, A::TBASIS_ELEM, B::TBASIS_ELEM)

    out = TBASIS_ELEM()

    for (alpha,c) in A
        for (beta,d) in B

            P = tbasis_mul_basis(Type, alpha, beta)

            for (gamma,e) in P
                out[gamma] = get(out, gamma, QQ(0)) + c*d*e
                if out[gamma] == 0
                    delete!(out, gamma)
                end
            end

        end
    end

    return out

end

#########

"""
    chebyshev(Type::Symbol, alpha::Vector{Int64})

Return the generalized Chebyshev polynomial indexed by `alpha`.

The input `alpha` is given in coordinates with respect to the basis of
fundamental weights. The polynomial is computed in variables z1,...,zn
over QQ.
"""
function chebyshev(Type::Symbol, alpha::Vector{Int64})

    n = length(alpha)
    R, z, cache = chebyshevcontext(Type,n)

    beta = pull(Type, alpha)
    key = Tuple(beta)

    if haskey(cache, key)
        return cache[key]
    end

    if all(x -> x == 0, beta)
        cache[key] = R(1)
        return cache[key]
    end

    for j in 1:n
        if beta[j] == 1 && all(beta[i] == 0 for i in 1:n if i != j)
            cache[key] = z[j]
            return cache[key]
        end
    end

    index = findfirst(x -> x > 0, beta)
    gamma = zeros(Int64, n)
    gamma[index] = 1
    base = beta .- gamma

    # use structure-constant cache
    M = cheb_mult(Type, base, gamma)  # Dict{Tuple{Int64},QQFieldElem}
    coeff_beta = get(M, Tuple(beta), QQ(0))
    @assert coeff_beta != 0

    # numerator: multiply T_gamma * T_base in z-variables
    numerator = chebyshev(Type, gamma) * chebyshev(Type, base)

    for (k, c) in M
        if k != Tuple(beta)
            numerator -= c * chebyshev(Type, collect(k))
        end
    end

    result = numerator / coeff_beta
    cache[key] = result

    return result
end

#########

function tbasis_to_z(Type::Symbol, F::TBASIS_ELEM)

    n = length(first(keys(F)))

    Rz, z, cache = chebyshevcontext(Type,n)

    out_z = Rz(0)

    for (alpha,c) in F
        out_z += Rz(c) * chebyshev(Type, collect(alpha))
    end

    return out_z

end

#########

function coefficient_matrices_T(M::Matrix{TBASIS_ELEM})

    r, c = size(M)

    coeffs = Dict{Tuple{Vararg{Int64}}, Matrix{QQFieldElem}}()

    for i in 1:r
        for j in 1:c
            for (alpha, q) in M[i,j]

                if !haskey(coeffs, alpha)
                    coeffs[alpha] = zeros(QQ, r, c)
                end

                coeffs[alpha][i,j] += q

            end
        end
    end

    return coeffs

end

# ============================================================
# 4. Multiplicative Invariants
# ============================================================

# Context cache for Laurent polynomial rings
const FUNDINV_CONTEXT = Dict{Int,Any}()

function fundinvcontext(n::Int)

    if haskey(FUNDINV_CONTEXT, n)
        return FUNDINV_CONTEXT[n]
    end

    names = vcat(
        [Symbol("x$i") for i in 1:n],
        [Symbol("y$i") for i in 1:n]
    )

    R, vars = polynomial_ring(QQ, names)

    xA = vars[1:n]
    yA = vars[n+1:2*n]

    FUNDINV_CONTEXT[n] = (R,xA,yA)

    return FUNDINV_CONTEXT[n]

end

# Fundamental invariants
function fundamental_invariant(Type::Symbol, n::Int)
    R, xA, yA = fundinvcontext(n)

    if Type == :A
    L = [xA[1]; [xA[k]*yA[k-1] for k in 2:n]; yA[n]]
    return [(1//binomial(n+1,l)) * esp(L,l) for l in 1:n]

    elseif Type == :C
        L = [xA[1]+yA[1]; [xA[k]*yA[k-1] + yA[k]*xA[k-1] for k in 2:n]]
        return [(1//(binomial(n,l)*2^l)) * esp(L,l) for l in 1:n]

    elseif Type == :B
        L = [xA[1]+yA[1]; [xA[k]*yA[k-1] + yA[k]*xA[k-1] for k in 2:n-1]; xA[n]^2*yA[n-1] + yA[n]^2*xA[n-1]]
        orb = orbit(Type, vcat(zeros(Int64,n-1), [1]))
        orb_terms = [prod((v[i] >= 0 ? xA[i]^v[i] : yA[i]^(-v[i])) for i in 1:n) for v in orb]
        return vcat(
            [(1//(binomial(n,l)*2^l)) * esp(L,l) for l in 1:n-1],
            [(1//2^n) * sum(orb_terms)]
        )

    elseif Type == :D
        L = [xA[1]+yA[1]; [xA[k]*yA[k-1] + yA[k]*xA[k-1] for k in 2:n-2];
             xA[n]*xA[n-1]*yA[n-2] + yA[n]*yA[n-1]*xA[n-2];
             xA[n]*yA[n-1] + yA[n]*xA[n-1]]
        orb1 = orbit(Type, vcat(zeros(Int64,n-2), [1,0]))
        orb2 = orbit(Type, vcat(zeros(Int64,n-1), [1]))

        orb_terms1 = [prod(v[i] >= 0 ? xA[i]^v[i] : yA[i]^(-v[i]) for i in 1:n) for v in orb1]
        orb_terms2 = [prod(v[i] >= 0 ? xA[i]^v[i] : yA[i]^(-v[i]) for i in 1:n) for v in orb2]

        return vcat(
            [(1//(binomial(n,l)*2^l)) * esp(L,l) for l in 1:n-2],
            [(1//2^(n-1)) * sum(orb_terms1)],
            [(1//2^(n-1)) * sum(orb_terms2)]
        )

    elseif Type == :F && n == 4

        T1 = 1//24*(xA[1]^2*yA[2] + xA[1]*xA[2]*yA[3]^2 + xA[1]*xA[3]^2*yA[2]^2 + xA[1]*xA[3]^2*yA[2]*yA[4]^2 + xA[1]*xA[4]^2*yA[2] + xA[1]*xA[4]^2*yA[3]^2 + xA[1]*yA[2] + xA[1]*yA[4]^2 + xA[1] + xA[2]^2*yA[1]*yA[3]^2 + xA[2]*xA[4]^2*yA[1]*yA[3]^2 + xA[2]*xA[4]^2*yA[3]^2 + xA[2]*yA[1]^2 + xA[2]*yA[1]*yA[4]^2 + xA[2]*yA[1] + xA[2]*yA[3]^2 + xA[2]*yA[4]^2 + xA[3]^2*yA[1]*yA[2] + xA[3]^2*yA[1]*yA[4]^2 + xA[3]^2*yA[2]*yA[4]^2 + xA[3]^2*yA[2] + xA[4]^2*yA[1] + xA[4]^2*yA[2] + yA[1])

        T2 = 1//96*(xA[1]^3*xA[3]^2*yA[2]^3 + xA[1]^3*xA[3]^2*yA[2]^2*yA[4]^2 + xA[1]^3*xA[4]^2*yA[2]^2 + xA[1]^3*xA[4]^2*yA[2]*yA[3]^2 + xA[1]^3*yA[2]^2 + xA[1]^3*yA[2]*yA[4]^2 + xA[1]^3*yA[2] + xA[1]^3*yA[3]^2 + xA[1]^2*xA[2]*xA[4]^2*yA[3]^4 + xA[1]^2*xA[2]*yA[3]^2*yA[4]^2 + xA[1]^2*xA[2]*yA[3]^2 + xA[1]^2*xA[3]^4*yA[2]^3*yA[4]^2 + xA[1]^2*xA[3]^2*xA[4]^2*yA[2]^3 + xA[1]^2*xA[3]^2*yA[2]^3 + xA[1]^2*xA[3]^2*yA[2]*yA[4]^4 + xA[1]^2*xA[3]^2*yA[2]*yA[4]^2 + xA[1]^2*xA[4]^4*yA[2]*yA[3]^2 + xA[1]^2*xA[4]^2*yA[2]*yA[3]^2 + xA[1]^2*xA[4]^2*yA[2] + xA[1]^2*yA[2]*yA[4]^2 + xA[1]*xA[2]^2*xA[4]^2*yA[3]^4 + xA[1]*xA[2]^2*yA[3]^4 + xA[1]*xA[2]^2*yA[3]^2*yA[4]^2 + xA[1]*xA[2]*xA[4]^4*yA[3]^4 + xA[1]*xA[2]*xA[4]^2*yA[3]^4 + xA[1]*xA[2]*xA[4]^2*yA[3]^2 + xA[1]*xA[2]*yA[3]^2*yA[4]^2 + xA[1]*xA[2]*yA[4]^4 + xA[1]*xA[2]*yA[4]^2 + xA[1]*xA[3]^4*yA[2]^3*yA[4]^2 + xA[1]*xA[3]^4*yA[2]^3 + xA[1]*xA[3]^4*yA[2]^2*yA[4]^4 + xA[1]*xA[3]^4*yA[2]^2*yA[4]^2 + xA[1]*xA[3]^2*xA[4]^2*yA[2]^3 + xA[1]*xA[3]^2*xA[4]^2*yA[2]^2 + xA[1]*xA[3]^2*yA[2]^2*yA[4]^2 + xA[1]*xA[3]^2*yA[2]*yA[4]^4 + xA[1]*xA[3]^2*yA[2] + xA[1]*xA[3]^2*yA[4]^4 + xA[1]*xA[4]^4*yA[2]^2 + xA[1]*xA[4]^4*yA[2]*yA[3]^2 + xA[1]*xA[4]^4*yA[3]^2 + xA[1]*xA[4]^2*yA[2]^2 + xA[1]*yA[3]^2 + xA[2]^3*xA[4]^2*yA[1]^2*yA[3]^4 + xA[2]^3*xA[4]^2*yA[1]*yA[3]^4 + xA[2]^3*yA[1]^3*yA[3]^2 + xA[2]^3*yA[1]^2*yA[3]^2*yA[4]^2 + xA[2]^3*yA[1]^2*yA[3]^2 + xA[2]^3*yA[1]*yA[3]^4 + xA[2]^3*yA[1]*yA[3]^2*yA[4]^2 + xA[2]^3*yA[3]^4 + xA[2]^2*xA[4]^4*yA[1]*yA[3]^4 + xA[2]^2*xA[4]^2*yA[1]^3*yA[3]^2 + xA[2]^2*xA[4]^2*yA[1]*yA[3]^4 + xA[2]^2*xA[4]^2*yA[1]*yA[3]^2 + xA[2]^2*yA[1]^3*yA[4]^2 + xA[2]^2*yA[1]^3 + xA[2]^2*yA[1]*yA[3]^2*yA[4]^2 + xA[2]^2*yA[1]*yA[4]^4 + xA[2]^2*yA[1]*yA[4]^2 + xA[2]*xA[3]^2*yA[1]^3*yA[4]^2 + xA[2]*xA[3]^2*yA[1]^2*yA[4]^4 + xA[2]*xA[3]^2*yA[1]^2*yA[4]^2 + xA[2]*xA[3]^2*yA[1]*yA[4]^4 + xA[2]*xA[4]^4*yA[1]^2*yA[3]^2 + xA[2]*xA[4]^4*yA[1]*yA[3]^2 + xA[2]*xA[4]^4*yA[3]^4 + xA[2]*xA[4]^2*yA[1]^3 + xA[2]*xA[4]^2*yA[1]^2*yA[3]^2 + xA[2]*xA[4]^2*yA[1]^2 + xA[2]*yA[1]^3 + xA[2]*yA[1]^2*yA[4]^2 + xA[2]*yA[1]*yA[3]^2 + xA[2]*yA[4]^4 + xA[2] + xA[3]^4*yA[1]^2*yA[2]*yA[4]^2 + xA[3]^4*yA[1]*yA[2]^2*yA[4]^2 + xA[3]^4*yA[1]*yA[2]^2 + xA[3]^4*yA[1]*yA[2]*yA[4]^4 + xA[3]^4*yA[1]*yA[2]*yA[4]^2 + xA[3]^4*yA[2]^3 + xA[3]^4*yA[2]*yA[4]^4 + xA[3]^2*xA[4]^2*yA[1]^2*yA[2] + xA[3]^2*xA[4]^2*yA[1]*yA[2]^2 + xA[3]^2*xA[4]^2*yA[1]*yA[2] + xA[3]^2*yA[1]^3 + xA[3]^2*yA[1]^2*yA[2] + xA[3]^2*yA[1]*yA[2]*yA[4]^2 + xA[3]^2*yA[1]*yA[4]^4 + xA[3]^2*yA[1] + xA[4]^4*yA[1]*yA[2] + xA[4]^4*yA[1]*yA[3]^2 + xA[4]^4*yA[2] + xA[4]^2*yA[1]*yA[2] + yA[2])

        T3 = 1//96*(xA[1]^2*xA[3]^2*yA[2]^2*yA[4] + xA[1]^2*xA[3]*xA[4]*yA[2]^2 + xA[1]^2*xA[3]*yA[2]^2 + xA[1]^2*xA[3]*yA[2]*yA[4]^2 + xA[1]^2*xA[3]*yA[2]*yA[4] + xA[1]^2*xA[4]^2*yA[2]*yA[3] + xA[1]^2*xA[4]*yA[2]*yA[3] + xA[1]^2*xA[4]*yA[2] + xA[1]^2*xA[4]*yA[3]^2 + xA[1]^2*yA[2]*yA[4] + xA[1]^2*yA[3]*yA[4] + xA[1]^2*yA[3] + xA[1]*xA[2]*xA[4]^2*yA[3]^3 + xA[1]*xA[2]*xA[4]*yA[3]^3 + xA[1]*xA[2]*xA[4]*yA[3]^2 + xA[1]*xA[2]*yA[3]^2*yA[4] + xA[1]*xA[2]*yA[3]*yA[4]^2 + xA[1]*xA[2]*yA[3]*yA[4] + xA[1]*xA[3]^3*yA[2]^2*yA[4]^2 + xA[1]*xA[3]^3*yA[2]^2*yA[4] + xA[1]*xA[3]^2*xA[4]*yA[2]^2 + xA[1]*xA[3]^2*yA[2]^2*yA[4] + xA[1]*xA[3]^2*yA[2]*yA[4]^3 + xA[1]*xA[3]^2*yA[2]*yA[4] + xA[1]*xA[3]*xA[4]^2*yA[2]^2 + xA[1]*xA[3]*xA[4]*yA[2]^2 + xA[1]*xA[3]*xA[4]*yA[2] + xA[1]*xA[3]*yA[2]*yA[4]^2 + xA[1]*xA[3]*yA[4]^3 + xA[1]*xA[3]*yA[4]^2 + xA[1]*xA[4]^3*yA[2]*yA[3] + xA[1]*xA[4]^3*yA[3]^2 + xA[1]*xA[4]^2*yA[2]*yA[3] + xA[1]*xA[4]^2*yA[3] + xA[1]*xA[4]*yA[3]^2 + xA[1]*yA[3]*yA[4] + xA[2]^2*xA[4]^2*yA[1]*yA[3]^3 + xA[2]^2*xA[4]*yA[1]^2*yA[3]^2 + xA[2]^2*xA[4]*yA[1]*yA[3]^3 + xA[2]^2*xA[4]*yA[1]*yA[3]^2 + xA[2]^2*xA[4]*yA[3]^3 + xA[2]^2*yA[1]^2*yA[3]*yA[4] + xA[2]^2*yA[1]^2*yA[3] + xA[2]^2*yA[1]*yA[3]^2*yA[4] + xA[2]^2*yA[1]*yA[3]*yA[4]^2 + xA[2]^2*yA[1]*yA[3]*yA[4] + xA[2]^2*yA[3]^3 + xA[2]^2*yA[3]^2*yA[4] + xA[2]*xA[3]*yA[1]^2*yA[4]^2 + xA[2]*xA[3]*yA[1]^2*yA[4] + xA[2]*xA[3]*yA[1]*yA[4]^3 + xA[2]*xA[3]*yA[1]*yA[4]^2 + xA[2]*xA[4]^3*yA[1]*yA[3]^2 + xA[2]*xA[4]^3*yA[3]^3 + xA[2]*xA[4]^2*yA[1]^2*yA[3] + xA[2]*xA[4]^2*yA[1]*yA[3] + xA[2]*xA[4]^2*yA[3]^3 + xA[2]*xA[4]*yA[1]^2*yA[3] + xA[2]*xA[4]*yA[1]^2 + xA[2]*xA[4]*yA[1]*yA[3]^2 + xA[2]*xA[4]*yA[3] + xA[2]*yA[1]^2*yA[4] + xA[2]*yA[1]*yA[3]*yA[4] + xA[2]*yA[3]*yA[4]^2 + xA[2]*yA[4]^3 + xA[2]*yA[4] + xA[3]^3*yA[1]*yA[2]*yA[4]^2 + xA[3]^3*yA[1]*yA[2]*yA[4] + xA[3]^3*yA[2]^2*yA[4] + xA[3]^3*yA[2]^2 + xA[3]^3*yA[2]*yA[4]^3 + xA[3]^3*yA[2]*yA[4]^2 + xA[3]^2*xA[4]*yA[1]*yA[2] + xA[3]^2*xA[4]*yA[2]^2 + xA[3]^2*yA[1]^2*yA[4] + xA[3]^2*yA[1]*yA[2]*yA[4] + xA[3]^2*yA[1]*yA[4]^3 + xA[3]^2*yA[1]*yA[4] + xA[3]^2*yA[4]^3 + xA[3]*xA[4]^2*yA[1]*yA[2] + xA[3]*xA[4]^2*yA[2] + xA[3]*xA[4]*yA[1]^2 + xA[3]*xA[4]*yA[1]*yA[2] + xA[3]*xA[4]*yA[1] + xA[3]*yA[1]^2 + xA[3]*yA[1]*yA[4]^2 + xA[3]*yA[2]*yA[4] + xA[3]*yA[4]^3 + xA[3] + xA[4]^3*yA[1]*yA[3] + xA[4]^3*yA[2] + xA[4]^3*yA[3]^2 + xA[4]^3*yA[3] + xA[4]^2*yA[1]*yA[3] + xA[4]*yA[2] + yA[3])

        T4 = 1//24*(xA[1]*xA[3]*yA[2]*yA[4] + xA[1]*xA[3]*yA[2] + xA[1]*xA[4]*yA[2] + xA[1]*xA[4]*yA[3] + xA[1]*yA[3] + xA[1]*yA[4] + xA[2]*xA[4]*yA[1]*yA[3] + xA[2]*xA[4]*yA[3]^2 + xA[2]*yA[1]*yA[3] + xA[2]*yA[1]*yA[4] + xA[2]*yA[3]*yA[4] + xA[2]*yA[3] + xA[3]^2*yA[2]*yA[4] + xA[3]*xA[4]*yA[2] + xA[3]*yA[1]*yA[4] + xA[3]*yA[1] + xA[3]*yA[2] + xA[3]*yA[4]^2 + xA[3]*yA[4] + xA[4]^2*yA[3] + xA[4]*yA[1] + xA[4]*yA[3] + xA[4] + yA[4])

    return [T1,T2,T3,T4]

    elseif Type == :E && 6 <= n <= 8
        invariants = []
        for i in 1:n
            ei = [j==i ? 1 : 0 for j in 1:n]
            orb = orbit(Type,ei)
            ord = orbitcardinality(Type,ei)
            s = (1//ord) * sum(
                prod(v[k] >= 0 ? xA[k]^v[k] : yA[k]^(-v[k]) for k in 1:n)
                for v in orb
            )
            push!(invariants,s)
        end
        return invariants

    elseif Type == :G && n == 2
    T1 = (1//6) * (
        xA[1] +
        xA[2]*yA[1] +
        xA[1]^2*yA[2] +
        xA[2]*yA[1]^2 +
        xA[1]*yA[2] +
        yA[1]
    )
    T2 = (1//6) * (
        xA[2] +
        xA[2]^2*yA[1]^3 +
        xA[1]^3*yA[2] +
        xA[1]^3*yA[2]^2 +
        xA[2]*yA[1]^3 +
        yA[2]
    )
    return [T1,T2]

    else
        error("Root system must be simple type A, B, C, D, E, F, G")
    end

end

#############

function multiplicative_action(Type::Symbol, s::Matrix{Int}, f::QQMPolyRingElem)
    R = parent(f)
    vars = gens(R) 
    n = div(length(vars),2)
    out = R(0)
    mons = collect(monomials(f))
    coeffs = collect(coefficients(f))

    for (m,c) in zip(mons, coeffs)
        e = collect(exponents(m))[1]
        alpha = [Int(e[i] - e[n+i]) for i in 1:n]  # exponent vector in fundamental weights
        # compute the image of exponents under the integer matrix s
        alpha_image = s * alpha
        # reconstruct the Laurent monomial
        m_image = prod(alpha_image[i] >= 0 ? vars[i]^alpha_image[i] : vars[n+i]^(-alpha_image[i]) for i in 1:n)
        out += QQ(c) * m_image
    end

    return out
end

function is_multiplicative_invariant(Type::Symbol, f::QQMPolyRingElem)
    R = parent(f)
    n = div(length(gens(R)), 2)
    W = zweylgroupgen(Type, n)

    for s in W
        f_image = multiplicative_action(Type, s, f)
        if f != f_image
            return false
        end
    end

    return true
end

function invariant_rewrite_T(Type::Symbol, f::QQMPolyRingElem)

    if !is_multiplicative_invariant(Type, f)
        error("Input polynomial is not a multiplicative invariant")
    end

    R = parent(f)
    vars = gens(R) 
    n = div(length(vars),2)

    out_T = TBASIS_ELEM()

    mons = collect(monomials(f))
    coeffs = collect(coefficients(f))

    for (m,c) in zip(mons, coeffs)

        e = collect(exponents(m))[1]

        alpha = [Int64(e[i] - e[n+i]) for i in 1:n]

        if all(a -> a >= 0, alpha)

            q = QQ(c) * QQ(orbitcardinality(Type, alpha))

            key = Tuple(pull(Type, alpha))
            out_T[key] = get(out_T, key, QQ(0)) + q

            if out_T[key] == 0
                delete!(out_T, key)
            end

        end

    end

    return out_T

end

function invariant_rewrite_z(Type::Symbol, f::QQMPolyRingElem)

    if !is_multiplicative_invariant(Type, f)
        error("Input polynomial is not a multiplicative invariant")
    end

    return tbasis_to_z(Type,invariant_rewrite_T(Type,f))

end

function invariant_rewrite(Type::Symbol, f::QQMPolyRingElem)

    if !is_multiplicative_invariant(Type, f)
        error("Input polynomial is not a multiplicative invariant")
    end

    out_T = invariant_rewrite_T(Type,f)
    out_z = invariant_rewrite_z(Type,f)

    return (out_T,out_z)

end

#############



# ============================================================
# 5. Moments
# ============================================================

function hermite_entry_T(Type::Symbol, n::Int, k::Int)
    
    out = tbasis_term(Type, vcat(k, zeros(Int,n-1)))
    out = tbasis_scale(-1, out)

    if Type in (:B,:C,:D,:G)
        if isodd(k)
            for j in 1:div(k-1,2)
                alpha = tbasis_term(Type, vcat(k - 2*j, zeros(Int,n-1)))
                factor = 4*binomial(k-2,j-1) - binomial(k,j)
                alpha = tbasis_scale(factor, alpha)
                out = tbasis_add(out, alpha)
            end
        else
            for j in 1:(div(k,2)-1)
                alpha = tbasis_term(Type, vcat(k - 2*j, zeros(Int,n-1)))
                factor = 4*binomial(k-2,j-1) - binomial(k,j)
                alpha = tbasis_scale(factor, alpha)
                out = tbasis_add(out, alpha)
            end
            alpha = tbasis_term(Type, zeros(Int,n))
            factor = div(4*binomial(k-2,div(k,2)-1) - binomial(k,div(k,2)), 2)
            alpha = tbasis_scale(factor, alpha)
            out = tbasis_add(out, alpha)
        end
    else
        error("Root system must be of simple type B,C,D,G")
    end

    return out
end

function hermite_matrix_T(Type::Symbol, n::Int)
    if Type in (:B, :C, :D)
        dim = n
    elseif Type == :G
        dim = 3
    else
        error("Root system must be B, C, D, or G2 for this Hermite matrix")
    end

    M = Matrix{TBASIS_ELEM}(undef, dim, dim)

    for i in 1:dim, j in 1:dim
        M[i,j] = hermite_entry_T(Type, n, i + j)
    end

    return M
end

#############

function moment_matrix_T(Type::Symbol, n::Int, d::Int)

    BT = chebyshevbasis_upto(Type,n,d)

    N = length(BT)

    MT = Matrix{TBASIS_ELEM}(undef,N,N)

    for i in 1:N
        for j in 1:N
            MT[i,j] = tbasis_mul_basis(Type, Tuple(BT[i]), Tuple(BT[j]))
        end
    end

    return MT

end

function localized_moment_matrix_T(Type::Symbol, n::Int, d::Int)

    if !(Type in (:B,:C,:D,:G))
        error("This function is implemented only for types B, C, D")
    end

    MT = moment_matrix_T(Type,n,d)

    N = size(MT,1)

    if Type in (:B,:C,:D)
        hdim = n
    elseif Type == :G 
        hdim = 3
    end

    LMT = Matrix{TBASIS_ELEM}(undef, hdim*N, hdim*N)

    for a in 1:hdim
        for b in 1:hdim
            for i in 1:N
                for j in 1:N

                    row = (a-1)*N + i
                    col = (b-1)*N + j

                    hT = hermite_entry_T(Type, n, a+b)
                    mT = MT[i,j]

                    LMT[row,col] = tbasis_mul(Type, hT, mT)

                end
            end
        end
    end

    return LMT

end

#############

using JuMP
using MosekTools
using SCS

function chromatic_sdp_bound(Type::Symbol, n::Int, d::Int, r::Int, solver)

    if Type in (:B, :C)
        D = n
    elseif Type == :D && iseven(n)
        D = n
    elseif Type == :G && n == 2
        D = 3
    else
        error("Implemented only for Bn, Cn, D2n, and G2")
    end

    if d < D
        error("Need d >= D")
    end

    M = moment_matrix_T(Type, n, d)
    LM = localized_moment_matrix_T(Type, n, d - D)

    A = coefficient_matrices_T(M)
    B = coefficient_matrices_T(LM)

    NA = size(first(values(A)), 1)
    NB = size(first(values(B)), 1)

    BT = chebyshevbasis_upto(Type, n, 2*d)
    all_indices = [Tuple(alpha) for alpha in BT]

    zero_key = Tuple(zeros(Int64,n))
    S = [Tuple(alpha) for alpha in chebyshevlevel(Type,n,r)]
    complement = setdiff(all_indices, union(S, Set([zero_key])))

    A0 = get(A, zero_key, zeros(QQ,NA,NA))
    B0 = get(B, zero_key, zeros(QQ,NB,NB))

    model = Model(solver.Optimizer)

    @variable(model, XA[1:NA,1:NA], PSD)
    @variable(model, XB[1:NB,1:NB], PSD)



    @objective(model, Max,
        - sum(Float64(A0[i,j]) * XA[i,j] for i in 1:NA, j in 1:NA)
        - sum(Float64(B0[i,j]) * XB[i,j] for i in 1:NB, j in 1:NB)
    )



    # Precompute the sums of A_alpha and B_alpha over alpha in S
    A_sum = zeros(QQ, NA, NA)
    B_sum = zeros(QQ, NB, NB)

    for alpha in S
        A_sum .+= get(A, alpha, zeros(QQ, NA, NA))
        B_sum .+= get(B, alpha, zeros(QQ, NB, NB))
    end

    # Add the constraint using the summed matrices
    @constraint(model,
        sum(Float64(A_sum[i,j]) * XA[i,j] for i in 1:NA, j in 1:NA) +
        sum(Float64(B_sum[i,j]) * XB[i,j] for i in 1:NB, j in 1:NB)
        == 1
    )



    for alpha in S

        Aalpha = get(A, alpha, zeros(QQ,NA,NA))
        Balpha = get(B, alpha, zeros(QQ,NB,NB))

        expr =
            sum(Float64(Aalpha[i,j]) * XA[i,j] for i in 1:NA, j in 1:NA) +
            sum(Float64(Balpha[i,j]) * XB[i,j] for i in 1:NB, j in 1:NB)

        @constraint(model, expr >= 0)

    end



    for beta in complement

        Abeta = get(A, beta, zeros(QQ,NA,NA))
        Bbeta = get(B, beta, zeros(QQ,NB,NB))

        expr =
            sum(Float64(Abeta[i,j]) * XA[i,j] for i in 1:NA, j in 1:NA) +
            sum(Float64(Bbeta[i,j]) * XB[i,j] for i in 1:NB, j in 1:NB)

        @constraint(model, expr == 0)

    end

    optimize!(model)

    opt = objective_value(model)

    return 1 - 1/opt

end

#############

#using LinearAlgebra

function chromatic_sdp_data(Type::Symbol, n::Int, d::Int, r::Int; folder::String="/mnt/c/Users/Dr. Tobias Metzlaff/Documents/Mathematics/Projects/YAND/SDP/ChebyJev/Chromatic")

    # Compute moment matrices
    if Type in (:B,:C)
        D = n
    elseif Type == :D && iseven(n)
        D = n
    elseif Type == :G && n == 2
        D = 3
    else
        error("Implemented only for Bn, Cn, D2n, and G2")
    end
    if d < D
        error("Need d >= D")
    end

    M = ChebyJev.moment_matrix_T(Type,n,d)
    LM = ChebyJev.localized_moment_matrix_T(Type,n,d-D)

    A = ChebyJev.coefficient_matrices_T(M)
    B = ChebyJev.coefficient_matrices_T(LM)

    NA = size(first(values(A)),1)
    NB = size(first(values(B)),1)

    BT = ChebyJev.chebyshevbasis_upto(Type,n,2*d)
    all_indices = [Tuple(alpha) for alpha in BT]

    zero_key = Tuple(zeros(Int64,n))
    S = [Tuple(alpha) for alpha in ChebyJev.chebyshevlevel(Type,n,r)]
    complement = setdiff(all_indices, union(S, Set([zero_key])))

    A0 = get(A, zero_key, zeros(QQ,NA,NA))
    B0 = get(B, zero_key, zeros(QQ,NB,NB))

    # ----------------------
    # Step 3: SDPA file
    # ----------------------

    function write_sdpa_block_entries(io, mat_idx::Int, block_idx::Int, M)

        rows, cols = size(M)

        for i in 1:rows
            for j in i:cols
                val = Float64(M[i,j])
                if val != 0.0
                    println(io, "$mat_idx $block_idx $i $j $val")
                end
            end
        end

    end

    filename = joinpath(folder, "Chromatic_$(Type)$(n)_$(d)_$(r).dat-s")

    open(filename, "w") do io

        nconstr = 1 + length(S) + length(complement)
        nblocks = 2 + length(S)

        println(io, nconstr)
        println(io, nblocks)
        println(io, join(vcat([NA, NB], fill(1, length(S))), " "))

        bvec = zeros(Float64, nconstr)
        bvec[1] = 1.0
        println(io, join(bvec, " "))

        A0 = get(A, zero_key, zeros(QQ, NA, NA))
        B0 = get(B, zero_key, zeros(QQ, NB, NB))

        write_sdpa_block_entries(io, 0, 1, -A0)
        write_sdpa_block_entries(io, 0, 2, -B0)

        constr_idx = 1

        # Constraint 1:
        # sum_{alpha in S} <A_alpha,XA> + <B_alpha,XB> = 1
        for alpha in S
            write_sdpa_block_entries(io, constr_idx, 1, get(A, alpha, zeros(QQ, NA, NA)))
            write_sdpa_block_entries(io, constr_idx, 2, get(B, alpha, zeros(QQ, NB, NB)))
        end

        constr_idx += 1

        for (sidx, alpha) in enumerate(S)

            write_sdpa_block_entries(io, constr_idx, 1, get(A, alpha, zeros(QQ, NA, NA)))
            write_sdpa_block_entries(io, constr_idx, 2, get(B, alpha, zeros(QQ, NB, NB)))

            slack_block = 2 + sidx
            println(io, "$constr_idx $slack_block 1 1 -1.0")

            constr_idx += 1

        end

        for beta in complement

            write_sdpa_block_entries(io, constr_idx, 1, get(A, beta, zeros(QQ, NA, NA)))
            write_sdpa_block_entries(io, constr_idx, 2, get(B, beta, zeros(QQ, NB, NB)))

            constr_idx += 1

        end

    end

    println("SDP data saved to ", filename)
    return filename
end

#############

#############

#############

#############

#############

#############

#############



end