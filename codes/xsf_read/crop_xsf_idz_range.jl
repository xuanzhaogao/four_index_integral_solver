using Printf
include("/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/reader.jl")

const INPUT = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002_upper_half_to_bottom.xsf"
const OUTPUT = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002_cropped.xsf"
const IDZ_START = 51
const IDZ_END = 150

function write_xsf(path::AbstractString, structure, datagrid)
    open(path, "w") do io
        println(io, "CRYSTAL")
        println(io, "PRIMVEC")
        for r in 1:3
            @printf(io, " %15.8f %15.8f %15.8f\n", structure.primvec[r, 1], structure.primvec[r, 2], structure.primvec[r, 3])
        end

        println(io, "PRIMCOORD")
        nat = length(structure.species)
        println(io, "  $(nat) 1")
        for a in 1:nat
            @printf(io, "%-3s %15.8f %15.8f %15.8f\n", structure.species[a], structure.positions[a, 1], structure.positions[a, 2], structure.positions[a, 3])
        end

        println(io)
        println(io, "BEGIN_BLOCK_DATAGRID_3D")
        println(io, "3D_field")
        println(io, "BEGIN_DATAGRID_3D_UNKNOWN")

        nx, ny, nz = datagrid.nx, datagrid.ny, datagrid.nz
        println(io, @sprintf(" %5d %5d %5d", nx, ny, nz))
        @printf(io, " %15.8f %15.8f %15.8f\n", datagrid.origin[1], datagrid.origin[2], datagrid.origin[3])
        @printf(io, " %15.8f %15.8f %15.8f\n", datagrid.A[1], datagrid.A[2], datagrid.A[3])
        @printf(io, " %15.8f %15.8f %15.8f\n", datagrid.B[1], datagrid.B[2], datagrid.B[3])
        @printf(io, " %15.8f %15.8f %15.8f\n", datagrid.C[1], datagrid.C[2], datagrid.C[3])

        data = vec(datagrid.values)
        i = 1
        n = length(data)
        while i <= n
            j = min(i + 5, n)
            for t in i:j
                @printf(io, " % .8E", data[t])
            end
            println(io)
            i = j + 1
        end

        println(io, "END_DATAGRID_3D")
        println(io, "END_BLOCK_DATAGRID_3D")
    end
end

function crop_datagrid_idz(dg, kz1::Int, kz2::Int)
    nx, ny, nz = dg.nx, dg.ny, dg.nz
    (1 <= kz1 <= kz2 <= nz) || error("invalid idz range: $(kz1):$(kz2), nz=$(nz)")

    new_values = copy(dg.values[:, :, kz1:kz2])
    new_nz = size(new_values, 3)

    # Preserve physical coordinates of retained points by shifting origin by the old k-offset.
    At, Bt, Ct = true_cell_vectors(dg)
    shift = (kz1 - 1) / nz
    new_origin = (
        dg.origin[1] + shift * Ct[1],
        dg.origin[2] + shift * Ct[2],
        dg.origin[3] + shift * Ct[3],
    )

    # Keep x/y axes unchanged; scale C so spacing convention matches the cropped nz.
    scale_c = (new_nz - 1) / (nz - 1)
    new_C = (
        dg.C[1] * scale_c,
        dg.C[2] * scale_c,
        dg.C[3] * scale_c,
    )

    return (
        nx = nx,
        ny = ny,
        nz = new_nz,
        origin = SVector(new_origin...),
        A = dg.A,
        B = dg.B,
        C = SVector(new_C...),
        values = new_values,
    )
end

function main()
    structure, dg = read_xsf(INPUT)
    dg2 = crop_datagrid_idz(dg, IDZ_START, IDZ_END)
    write_xsf(OUTPUT, structure, dg2)
    println("Wrote: ", OUTPUT)
    println("Old grid: ", (dg.nx, dg.ny, dg.nz), " -> New grid: ", (dg2.nx, dg2.ny, dg2.nz))
    println("idz kept: ", IDZ_START, ":", IDZ_END)
end

main()
