using BoundaryIntegral
import BoundaryIntegral as BI

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1   = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2   = joinpath(REF_DIR, "graphene_00002.xsf")
const COQUI   = joinpath(REF_DIR, "_coqui_crpa_loc.out")

# Slab geometry
const LZ = 3.35
const Z_CENTER = LZ / 2
const EPS_IN = 2.4
const EPS_OUT = 1.0
const L = 90.0

# Solver tolerances (reused from the bilayer protocol)
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

# Density-channel bandwidth sweep (used by Task 7 main; Task 6 smoke uses only Sharp)
const BANDWIDTHS = [0.05, 0.1, 0.2, 0.4]

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_k323201.csv")

function _solve_kwargs()
    return (
        n_quad = N_QUAD,
        edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL,
        gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER,
        max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
end

function _record_kwargs(mode_label, bandwidth)
    return (
        mode_label = mode_label, bandwidth = bandwidth,
        eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
        source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
    )
end

function smoke()
    println("Loading sources ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    println("Running density Sharp solve ...")
    density_specs = density_target_specs(src)
    result = solve_screened_mode(
        src.vs1, density_specs,
        L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
        _solve_kwargs()...,
    )
    println("  residual = ", result.residual,
            "  n_interface_points = ", result.n_interface_points)
    for pr in result.pair_results
        println("  ", pr.pair, "  u_total_ev = ", pr.u_total_ev,
                "  (u_int_ev = ", pr.u_int_ev, ", u_scatter_ev = ", pr.u_scatter_ev, ")")
    end
    coqui = parse_coqui_loc(COQUI)
    println("CoQui reference:")
    for ch in (:onsite, :nn)
        println("  ", ch, "  U_ijkl = ", coqui[ch].U_ijkl)
    end
    return result
end

# Run smoke when invoked directly (Task 6).
smoke()
