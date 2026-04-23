"""
Diagnose the near-field correction matrix:
- Condition number of C (Jacobi conversion)
- Magnitude of F entries
- Column-by-column accuracy of A_near vs A_direct
"""

using LinearAlgebra, FastGaussQuadrature, HCubature, Printf

beta = 0.6703

function jacobi_poly_vals!(v, beta, t)
    n = length(v); n == 0 && return v
    v[1] = 1.0; n == 1 && return v
    v[2] = ((beta+2)*t - beta)/2
    for m in 1:(n-2)
        abm = 2.0*m+beta; mp1 = m+1.0; mbp1 = m+beta+1.0
        a_c = (abm+1)*(abm+2)/(2*mp1*mbp1)
        b_c = -(beta^2)*(abm+1)/(2*mp1*mbp1*abm)
        c_c = m*(m+beta)*(abm+2)/(mp1*mbp1*abm)
        v[m+2] = (a_c*t + b_c)*v[m+1] - c_c*v[m]
    end
    return v
end

jacobi_norms_sq(beta, n) = [2.0^(beta+1)/(2.0*m+beta+1) for m in 0:(n-1)]

L_panel = 0.05
mid   = (L_panel/2, 0.0)
half  = (L_panel/2, 0.0)
Lhalf = L_panel/2

function kernel_DT(src, trg, nrm)
    d = (trg[1]-src[1], trg[2]-src[2])
    r2 = d[1]^2 + d[2]^2
    return -(d[1]*nrm[1]+d[2]*nrm[2])/(2π*r2)
end

# ── Table 1: condition number of C ─────────────────────────────────────────
println("n_quad  cond(C)")
for nq in [4, 8, 12, 16, 20, 24]
    gj_xs_r, gj_ws_r = gaussjacobi(nq, 0.0, beta)
    gj_xs = Float64.(gj_xs_r); gj_ws = Float64.(gj_ws_r)
    h = jacobi_norms_sq(beta, nq)
    P_mat = [jacobi_poly_vals!(zeros(nq), beta, gj_xs[j])[m] for j in 1:nq, m in 1:nq]
    C = (P_mat' .* reshape(gj_ws, 1, nq)) ./ h
    @printf("%-7d %.3e\n", nq, cond(C))
end

# ── Table 2: column accuracy of A_near for target near corner ──────────────
println()
println("Column accuracy of A_near vs A_direct (n_quad=12, target at varying d)")
println("d          max|A_near - A_near_ref|  max|A_direct|  max|delta=A_near-A_direct|")
nq = 12
gj_xs_r, gj_ws_r = gaussjacobi(nq, 0.0, beta)
gj_xs = Float64.(gj_xs_r); gj_ws = Float64.(gj_ws_r)
h = jacobi_norms_sq(beta, nq)
P_mat = [jacobi_poly_vals!(zeros(nq), beta, gj_xs[j])[m] for j in 1:nq, m in 1:nq]
C = (P_mat' .* reshape(gj_ws, 1, nq)) ./ h
pv2 = zeros(nq)

for d_trg in [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005]
    xk = (0.0, d_trg); nk = (1.0, 0.0)

    # F[m] via hcub (rtol=1e-10)
    F = zeros(nq)
    for m in 1:nq
        F[m], _ = hquadrature(
            t -> begin
                s = (mid[1]+t*half[1], mid[2]+t*half[2])
                K = kernel_DT(s, xk, nk)
                jacobi_poly_vals!(pv2, beta, t)
                K*(1+t)^beta*pv2[m]*Lhalf
            end, -1.0, 1.0; rtol=1e-10, atol=1e-14)
    end
    A_near = vec(F' * C)   # length nq

    # F_ref via hcub (rtol=1e-13)
    F_ref = zeros(nq)
    for m in 1:nq
        F_ref[m], _ = hquadrature(
            t -> begin
                s = (mid[1]+t*half[1], mid[2]+t*half[2])
                K = kernel_DT(s, xk, nk)
                jacobi_poly_vals!(pv2, beta, t)
                K*(1+t)^beta*pv2[m]*Lhalf
            end, -1.0, 1.0; rtol=1e-13, atol=1e-15)
    end
    A_near_ref = vec(F_ref' * C)

    # A_direct
    A_direct = [kernel_DT((mid[1]+gj_xs[j]*half[1], mid[2]+gj_xs[j]*half[2]), xk, nk) *
                gj_ws[j]*Lhalf for j in 1:nq]

    delta = A_near .- A_direct
    @printf("%.2e  %-22.4e  %-14.4e  %-14.4e\n",
        d_trg, maximum(abs.(A_near .- A_near_ref)), maximum(abs.(A_direct)), maximum(abs.(delta)))
end

# ── Table 3: how does correction magnitude vs correction error scale with nq ─
println()
println("n_quad  d=0.001:  max|delta|    max|A_near err|  ratio")
d_trg = 0.001; xk = (0.0, d_trg); nk = (1.0, 0.0)
for nq in [4, 8, 12, 16, 20]
    gj_xs_r, gj_ws_r = gaussjacobi(nq, 0.0, beta)
    gj_xs = Float64.(gj_xs_r); gj_ws = Float64.(gj_ws_r)
    h = jacobi_norms_sq(beta, nq)
    P_mat = [jacobi_poly_vals!(zeros(nq), beta, gj_xs[j])[m] for j in 1:nq, m in 1:nq]
    C = (P_mat' .* reshape(gj_ws, 1, nq)) ./ h
    pv2 = zeros(nq)

    F = zeros(nq)
    for m in 1:nq
        F[m], _ = hquadrature(
            t -> begin
                s = (mid[1]+t*half[1], mid[2]+t*half[2])
                K = kernel_DT(s, xk, nk); jacobi_poly_vals!(pv2, beta, t)
                K*(1+t)^beta*pv2[m]*Lhalf
            end, -1.0, 1.0; rtol=1e-10, atol=1e-14)
    end
    F_ref = zeros(nq)
    for m in 1:nq
        F_ref[m], _ = hquadrature(
            t -> begin
                s = (mid[1]+t*half[1], mid[2]+t*half[2])
                K = kernel_DT(s, xk, nk); jacobi_poly_vals!(pv2, beta, t)
                K*(1+t)^beta*pv2[m]*Lhalf
            end, -1.0, 1.0; rtol=1e-13, atol=1e-15)
    end
    A_near = vec(F' * C); A_near_ref = vec(F_ref' * C)
    A_direct = [kernel_DT((mid[1]+gj_xs[j]*half[1], mid[2]+gj_xs[j]*half[2]), xk, nk)*gj_ws[j]*Lhalf for j in 1:nq]
    delta = A_near .- A_direct
    err = maximum(abs.(A_near .- A_near_ref))
    @printf("%-7d %-12.4e  %-14.4e  %.2e\n", nq, maximum(abs.(delta)), err, err/max(maximum(abs.(delta)),1e-20))
end
