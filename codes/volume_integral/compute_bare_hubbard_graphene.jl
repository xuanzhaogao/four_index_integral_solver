using BoundaryIntegral, LinearAlgebra, FBCPoisson
import BoundaryIntegral as BI
using CSV, DataFrames

file = joinpath(@__DIR__, "data/hubbard_graphene.csv")
df = CSV.write(file, DataFrame(pair = [], N_FFT = [], U_raw = [], U_ev = []))

# ─── Input files ──────────────────────────────────────────────────────────────
graphene_1 = joinpath(@__DIR__, "../../density_data/graphene_00001_5x5x1_shifted.xsf")
graphene_2 = joinpath(@__DIR__, "../../density_data/graphene_00002_5x5x1_shifted.xsf")

# ─── Load and square wavefunctions → orbital densities ────────────────────────
structure_1, datagrid_1 = BI.read_xsf(graphene_1)
datagrid_1.values .*= datagrid_1.values   # psi -> |psi|^2
structure_2, datagrid_2 = BI.read_xsf(graphene_2)
datagrid_2.values .*= datagrid_2.values

# ─── Geometry ─────────────────────────────────────────────────────────────────
# The z-shift centres graphene at z ≈ 0 (cell C-vector = 15.92 Å, atoms at z ≈ 8 Å)
z_shift  = -7.920155482424242   # Å
a1_prim  = (2.465, 0.0, 0.0)    # primitive lattice vector a₁

# ─── VolumeSource objects ─────────────────────────────────────────────────────
# Orbital 1 at origin
vs1    = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, z_shift), tol = 1e-6)
# Orbital 2 at origin (nearest-neighbor A-B, d = 1.42 Å)
vs2    = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, z_shift), tol = 1e-6)
# Orbital 1 shifted by a₁ (next-nearest A-A, d = 2.465 Å)
vs1_a1 = BoundaryIntegral.VolumeSource(datagrid_1, shift = (a1_prim[1], a1_prim[2], z_shift), tol = 1e-6)
# Orbital 2 shifted by a₁ (A-B 2nd shell, d = √(2.465²+1.42²) ≈ 2.845 Å)
vs2_a1 = BoundaryIntegral.VolumeSource(datagrid_2, shift = (a1_prim[1], a1_prim[2], z_shift), tol = 1e-6)

# ─── Normalization ─────────────────────────────────────────────────────────────
# XSF stores psi with ∫|psi|² dV ≈ V_prim ≈ 84.26 Å³ (≠ 1).
# Divide by N_m * N_n to get the interaction between *normalized* orbital densities.
N1 = sum(vs1.density .* vs1.weights)
N2 = sum(vs2.density .* vs2.weights)
println("Normalization:  N1 = $N1 Å³,  N2 = $N2 Å³")

# ─── Physical conversion factor ───────────────────────────────────────────────
# lfbc3d kernel = 1/(4π|r|), physical Coulomb = e²/(4πε₀|r|)
# → U_phys [eV] = U_raw [Å⁻¹] × 4π × 14.3996 [eV·Å] / (N_m × N_n)
const e2_4pieps0 = 14.3996  # e²/(4πε₀) in eV·Å

function to_eV(raw, Na, Nb)
    return raw * 4π * e2_4pieps0 / (Na * Nb)
end

# ─── Compute U values ─────────────────────────────────────────────────────────
nffts = 32:32:256
for nfft in nffts
    tol  = 1e-6

    # U_00 : on-site (A-A, d = 0)
    pot00   = lfbc3d(nfft, vs1.positions, vs1.density .* vs1.weights, vs1.positions,    tol, 1)
    U_00_raw = sum(pot00 .* vs1.density .* vs1.weights)

    # U_01 : nearest-neighbor (A-B, d = 1.42 Å)
    pot01   = lfbc3d(nfft, vs1.positions, vs1.density .* vs1.weights, vs2.positions,    tol, 1)
    U_01_raw = sum(pot01 .* vs2.density .* vs2.weights)

    # U_02 : next-nearest-neighbor (A-A, d = 2.465 Å)
    pot02   = lfbc3d(nfft, vs1.positions, vs1.density .* vs1.weights, vs1_a1.positions, tol, 1)
    U_02_raw = sum(pot02 .* vs1_a1.density .* vs1_a1.weights)

    # U_03 : A-B 2nd shell (d = √(2.465²+1.42²) ≈ 2.845 Å)
    pot03   = lfbc3d(nfft, vs1.positions, vs1.density .* vs1.weights, vs2_a1.positions, tol, 1)
    U_03_raw = sum(pot03 .* vs2_a1.density .* vs2_a1.weights)

    CSV.write(file, DataFrame(pair = ["U_00", "U_01", "U_02", "U_03"], N_FFT = nfft, U_raw = [U_00_raw, U_01_raw, U_02_raw, U_03_raw], U_ev = [to_eV(U_00_raw, N1, N1), to_eV(U_01_raw, N1, N2), to_eV(U_02_raw, N1, N1), to_eV(U_03_raw, N1, N2)]), append = true)
end

# ─── Report ───────────────────────────────────────────────────────────────────
# d01 = 1.42
# d02 = 2.465
# d03 = sqrt(2.465^2 + 1.42^2)

# println("\n=== Bare Coulomb Hubbard parameters for graphene π orbitals ===")
# println("  U_00 = $(round(to_eV(U_00_raw, N1, N1), digits=4)) eV   (on-site A-A,           d = 0 Å)")
# println("  U_01 = $(round(to_eV(U_01_raw, N1, N2), digits=4)) eV   (nearest-neighbor A-B,  d = $d01 Å)")
# println("  U_02 = $(round(to_eV(U_02_raw, N1, N1), digits=4)) eV   (next-nearest A-A,      d = $d02 Å)")
# println("  U_03 = $(round(to_eV(U_03_raw, N1, N2), digits=4)) eV   (A-B 2nd shell,         d = $(round(d03, digits=3)) Å)")
# println("\n  Raw values (Å⁻¹):")
# println("    U_00_raw = $U_00_raw")
# println("    U_01_raw = $U_01_raw")
# println("    U_02_raw = $U_02_raw")
# println("    U_03_raw = $U_03_raw")
