using Printf
include("/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/reader.jl")

const INPUT = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002.xsf"
const OUTPUT = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002_upper_half_to_bottom.xsf"

"""
Shift the 3D DATAGRID by half along k (C direction):
new[:, :, 1:nz/2] = old[:, :, nz/2+1:nz]
new[:, :, nz/2+1:nz] = old[:, :, 1:nz/2]
"""
function shift_upper_half_to_bottom(values::Array{Float64, 3})
    nx, ny, nz = size(values)
    iseven(nz) || error("nz must be even to split exactly in half; got nz = $nz")
    return circshift(values, (0, 0, nz ÷ 2))
end

"""
Shift all atoms by -0.5 of the third lattice direction with periodic wrap in fractional coordinates.
"""
function shift_atoms_upper_half_to_bottom(structure)
    R = structure.primvec          # rows are lattice vectors a,b,c
    RT = transpose(R)              # cart = RT * frac
    nat = size(structure.positions, 1)
    new_positions = similar(structure.positions)
    for a in 1:nat
        r = structure.positions[a, :]
        frac = RT \ r
        frac[3] = mod(frac[3] - 0.5, 1.0)
        new_r = RT * frac
        new_positions[a, 1] = new_r[1]
        new_positions[a, 2] = new_r[2]
        new_positions[a, 3] = new_r[3]
    end
    return XSFStructure(structure.primvec, structure.species, new_positions)
end

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
        n = length(data)
        i = 1
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

function main()
    structure, dg = read_xsf(INPUT)
    shifted = shift_upper_half_to_bottom(dg.values)
    shifted_structure = shift_atoms_upper_half_to_bottom(structure)
    dg2 = (nx = dg.nx, ny = dg.ny, nz = dg.nz, origin = dg.origin, A = dg.A, B = dg.B, C = dg.C, values = shifted)
    write_xsf(OUTPUT, shifted_structure, dg2)
    println("Wrote: ", OUTPUT)
    println("Shape: ", size(shifted), ", shift along k by ", dg.nz ÷ 2)
    println("First atom old/new z: ", structure.positions[1, 3], " -> ", shifted_structure.positions[1, 3])
end

main()
