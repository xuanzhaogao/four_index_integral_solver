using LinearAlgebra

# ------------------------------------------------------------
# Gram matrix G[i,j] = 1 / (gammas[i] + gammas[j] + a0 + 1).
# ------------------------------------------------------------
function gram_matrix(a0::Real, gammas::AbstractVector{T}) where {T<:Real}
    if a0 <= -1
        error("Need a0 > -1.")
    end
    N = length(gammas)
    one_T = one(T)
    G = Matrix{T}(undef, N, N)
    for i in 1:N, j in 1:N
        e = gammas[i] + gammas[j] + T(a0)
        if e <= -one_T
            error("Non-integrable Gram entry at ($i,$j).")
        end
        G[i, j] = one_T / (e + one_T)
    end
    return G
end

# ------------------------------------------------------------
# Orthonormal generalized power basis Q_n = sum_k C[n,k] x^gammas[k]
# under weight x^a0 on [0,1], via Cholesky of the Gram matrix.
# ------------------------------------------------------------
function generalized_power_basis(
    a0::Real,
    gammas::AbstractVector{<:Real};
    rtol::Real = 1e-12,
    use_bigfloat::Bool = true,
)
    if use_bigfloat
        Gbig = gram_matrix(BigFloat(a0), collect(BigFloat, gammas))
        κ = cond(Float64.(Gbig))
        Lbig = Matrix(cholesky(Symmetric(Gbig)).L)
        Cbig = Lbig \ Matrix{BigFloat}(I, size(Gbig)...)
        return Cbig, Gbig, κ
    else
        G = gram_matrix(a0, collect(Float64, gammas))
        κ = cond(G)
        if κ > 1 / rtol
            @warn "Gram matrix ill-conditioned: cond(G) = $κ"
        end
        L = Matrix(cholesky(Symmetric(G)).L)
        C = L \ I(size(G, 1))
        return C, G, κ
    end
end


# ------------------------------------------------------------
# Evaluate Q_k(x) = sum_l C[k,l] x^gammas[l].
#
# Input:
#   x      : vector of nodes in (0,1)
#   gammas : exponents gamma_l
#   C      : coefficient matrix, rows are basis functions
#
# Output:
#   Qvals[k,j] = Q_k(x_j)
# ------------------------------------------------------------
function eval_Q_matrix(x::AbstractVector, gammas::AbstractVector, C::AbstractMatrix)
    nbasis = size(C, 1)
    npow   = length(gammas)
    nnode  = length(x)

    T = promote_type(eltype(x), eltype(gammas), eltype(C))
    Powers = Matrix{T}(undef, npow, nnode)

    for l in 1:npow
        γ = gammas[l]
        for j in 1:nnode
            Powers[l, j] = x[j]^γ
        end
    end

    return C * Powers
end


# ------------------------------------------------------------
# Evaluate derivative Q_k'(x).
#
# Q_k'(x) = sum_l C[k,l] gamma_l x^(gamma_l - 1).
#
# The gamma=0 term contributes zero.
# ------------------------------------------------------------
function eval_dQ_matrix(x::AbstractVector, gammas::AbstractVector, C::AbstractMatrix)
    nbasis = size(C, 1)
    npow   = length(gammas)
    nnode  = length(x)

    T = promote_type(eltype(x), eltype(gammas), eltype(C))
    dPowers = Matrix{T}(undef, npow, nnode)
    one_T = one(T)
    zero_T = zero(T)

    for l in 1:npow
        γ = gammas[l]

        if iszero(γ)
            for j in 1:nnode
                dPowers[l, j] = zero_T
            end
        else
            for j in 1:nnode
                dPowers[l, j] = γ * x[j]^(γ - one_T)
            end
        end
    end

    return C * dPowers
end


# ------------------------------------------------------------
# Logistic map to keep nodes inside (0,1).
# ------------------------------------------------------------
σ(y) = one(y) / (one(y) + exp(-y))

function logistic_vec(y::AbstractVector)
    return [σ(yi) for yi in y]
end

function logit(x::Real)
    return log(x / (one(x) - x))
end


# ------------------------------------------------------------
# Construct the exact moments of the orthonormal Q basis.
#
# If Q_0 = sqrt(a0+1), then
#
#   μ_0 = 1 / sqrt(a0+1),
#   μ_k = 0 for k >= 1.
#
# More generally, we compute μ_k from coefficients:
#
#   μ_k = sum_l C[k,l] / (a0 + gammas[l] + 1).
# ------------------------------------------------------------
function moments_Q(a0::Real, gammas::AbstractVector, C::AbstractMatrix, M::Integer)
    T = promote_type(typeof(a0), eltype(gammas), eltype(C))
    μ = zeros(T, M)
    one_T = one(T)

    for k in 1:M
        s = zero(T)
        for l in 1:length(gammas)
            s += C[k, l] / (T(a0) + gammas[l] + one_T)
        end
        μ[k] = s
    end

    return μ
end


# ------------------------------------------------------------
# Generalized Gaussian quadrature by Newton iteration.
#
# We solve:
#
#   sum_j w_j Q_k(x_j) = μ_k,   k=0,...,2n-1.
#
# Unknowns:
#   y_j, where x_j = logistic(y_j)
#   w_j
#
# The unknown vector is z = [y_1,...,y_n, w_1,...,w_n].
#
# Input:
#   a0     : weight exponent x^a0
#   gammas : exponents for generalized powers
#   C      : coefficient matrix for Q basis; must have at least 2n rows
#   n      : number of quadrature nodes
#
# Output:
#   x, w
# ------------------------------------------------------------
function generalized_gauss_quadrature(
    a0::Real,
    gammas::AbstractVector,
    C::AbstractMatrix,
    n::Integer;
    maxit::Integer = 30,
    tol::Real = 1e-13,
    verbose::Bool = true,
    x_init::Union{Nothing,AbstractVector} = nothing,
    w_init::Union{Nothing,AbstractVector} = nothing,
)
    M = 2n

    if size(C, 1) < M
        error("Need at least 2n orthonormal basis functions in C.")
    end

    # Use only Q_0,...,Q_{2n-1}.
    Cuse = C[1:M, :]

    # Target moments.
    μ = moments_Q(a0, gammas, Cuse, M)

    T = eltype(μ)

    if x_init === nothing
        x0 = T[T(0.5) * (one(T) - cos(T(2j - 1) * T(π) / T(2n))) for j in 1:n]
    else
        x0 = collect(T, x_init)
    end

    if w_init === nothing
        Qn = eval_Q_matrix(x0, gammas, Cuse[1:n, :])
        μn = μ[1:n]
        w0 = Qn \ μn
    else
        w0 = collect(T, w_init)
    end

    # Convert nodes to unconstrained variables y.
    y = T[logit(x0[j]) for j in 1:n]
    w = copy(w0)

    for it in 1:maxit
        x = logistic_vec(y)

        Q  = eval_Q_matrix(x, gammas, Cuse)
        dQ = eval_dQ_matrix(x, gammas, Cuse)

        # Residual F_k = sum_j w_j Q_k(x_j) - μ_k.
        F = Q * w - μ

        res = norm(F, Inf)

        if verbose
            println("Newton iter $it: residual = $res")
        end

        if res < tol
            return x, w
        end

        J = zeros(T, M, 2n)

        for k in 1:M
            for j in 1:n
                J[k, j] = w[j] * dQ[k, j] * x[j] * (one(T) - x[j])
                J[k, n + j] = Q[k, j]
            end
        end

        δ = -J \ F

        λ = one(T)
        accepted = false

        for ls in 1:60
            ytrial = y .+ λ .* δ[1:n]
            wtrial = w .+ λ .* δ[n+1:2n]

            xtrial = logistic_vec(ytrial)
            Qtrial = eval_Q_matrix(xtrial, gammas, Cuse)
            Ftrial = Qtrial * wtrial - μ

            if norm(Ftrial, Inf) < res
                y = ytrial
                w = wtrial
                accepted = true
                break
            end

            λ *= T(0.5)
        end

        if !accepted
            error("Newton line search failed.")
        end
    end

    error("Newton did not converge within maxit iterations.")
end


# ------------------------------------------------------------
# Continuation / homotopy solver.
#
# γ(λ) = (1-λ) γ_poly + λ γ_target,   γ_poly = [0,1,2,...,2m-1].
#
# At λ=0 the basis is the standard orthonormal polynomial basis on
# [0,1] under weight x^a0, and the m-point rule is the standard
# Gauss-Jacobi rule (mapped from [-1,1]).  We march λ from 0 → 1,
# rebuilding the basis and running Newton at each step warm-started
# from the previous solution.
#
# Returns (x, w) for the target gammas at λ=1.
# ------------------------------------------------------------
function generalized_gauss_continuation(
    a0::Real,
    gammas_target::AbstractVector,
    m::Integer;
    nsteps::Integer = 20,
    T = BigFloat,
    verbose::Bool = false,
    inner_maxit::Integer = 200,
    inner_tol::Real = T(1e-40),
)
    γ_target = collect(T, gammas_target[1:2m])
    γ_poly   = T[T(k - 1) for k in 1:2m]   # [0, 1, 2, ..., 2m-1]

    # Initial (x, w): standard Gauss-Jacobi for x^a0 on [0,1].
    a0T = T(a0)
    tk_f, Wk_f = gaussjacobi(m, 0.0, Float64(a0))
    x = T[T(0.5) * (one(T) + T(tk_f[k])) for k in 1:m]
    w = T[T(0.5)^(a0T + one(T)) * T(Wk_f[k]) for k in 1:m]

    λ_vals = range(zero(T), one(T); length = nsteps + 1)

    for (step, λ) in enumerate(λ_vals)
        γ = (one(T) - λ) .* γ_poly .+ λ .* γ_target
        C, _, κ = generalized_power_basis(a0T, γ; use_bigfloat = (T == BigFloat))
        if T != BigFloat
            C = T.(C)
        end

        if step == 1
            # λ=0 is already exactly satisfied by GJ rule; skip Newton.
            if verbose
                println("step 1 (λ=0)  cond(G)=", κ, "  using GJ exactly.")
            end
            continue
        end

        try
            x, w = generalized_gauss_quadrature(
                a0T, γ, C, m;
                x_init = x, w_init = w,
                maxit = inner_maxit, tol = inner_tol, verbose = false,
            )
        catch e
            error("Continuation failed at step $step (λ=$(Float64(λ))): $e")
        end

        if verbose
            println("step $step  λ=", Float64(λ),
                    "  cond(G)=", κ,
                    "  min(w)=", Float64(minimum(w)))
        end
    end

    return x, w
end