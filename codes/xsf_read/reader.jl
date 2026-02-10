using LinearAlgebra, StaticArrays
using Statistics
using GLMakie

struct XSFStructure
    primvec::Matrix{Float64}          # 3x3, rows are lattice vectors (Å or bohr depending on source)
    species::Vector{String}           # length nat
    positions::Matrix{Float64}        # nat x 3 Cartesian coordinates
end

struct XSFDatagrid3D
    nx::Int
    ny::Int
    nz::Int
    origin::SVector{3,Float64}
    A::SVector{3,Float64}
    B::SVector{3,Float64}
    C::SVector{3,Float64}
    values::Array{Float64,3}          # size (nx, ny, nz), i-fastest
end

S3(x::AbstractVector{<:Real}) = StaticArrays.SVector{3,Float64}(Float64(x[1]), Float64(x[2]), Float64(x[3]))

function _read_vec3(io)::Vector{Float64}
    while !eof(io)
        s = strip(readline(io))
        isempty(s) && continue
        return parse.(Float64, split(replace(s, "D" => "E")))
    end
    error("Unexpected EOF while reading vec3")
end

function read_xsf(path::AbstractString)
    primvec = nothing
    species = String[]
    positions = Matrix{Float64}(undef, 0, 3)

    datagrid = nothing

    open(path, "r") do io
        lines = readlines(io)
        nlines = length(lines)

        i = 1
        while i <= nlines
            s = strip(lines[i])

            if s == "PRIMVEC"
                v1 = parse.(Float64, split(strip(lines[i+1])))
                v2 = parse.(Float64, split(strip(lines[i+2])))
                v3 = parse.(Float64, split(strip(lines[i+3])))
                primvec = vcat(v1', v2', v3')  # 3x3 rows = vectors
                i += 4
                continue
            end

            if s == "PRIMCOORD"
                nat = parse(Int, split(strip(lines[i+1]))[1])
                species = Vector{String}(undef, nat)
                positions = Matrix{Float64}(undef, nat, 3)
                for a = 1:nat
                    parts = split(strip(lines[i+1+a]))
                    species[a] = parts[1]
                    positions[a, :] = parse.(Float64, parts[2:4])
                end
                i += 2 + nat
                continue
            end

            if startswith(s, "BEGIN_DATAGRID_3D")
                # header
                dims = parse.(Int, split(strip(lines[i+1])))
                nx, ny, nz = dims
                origin = S3(parse.(Float64, split(strip(lines[i+2]))))
                A = S3(parse.(Float64, split(strip(lines[i+3]))))
                B = S3(parse.(Float64, split(strip(lines[i+4]))))
                C = S3(parse.(Float64, split(strip(lines[i+5]))))

                n = nx * ny * nz
                data = Vector{Float64}(undef, n)
                p = 1
                j = i + 6
                while j <= nlines
                    t = strip(lines[j])
                    if startswith(t, "END_DATAGRID_3D")
                        break
                    end
                    if isempty(t)
                        j += 1
                        continue
                    end
                    parts = split(replace(t, "D" => "E"))   # handle Fortran exponent
                    vals = parse.(Float64, parts)
                    data[p:p+length(vals)-1] .= vals
                    p += length(vals)
                    j += 1
                    p > n && break
                end
                p == n+1 || error("DATAGRID expected $n values, got $(p-1)")

                values = reshape(data, (nx, ny, nz)) # i (x) fastest, then j (y), then k (z)
                datagrid = (nx=nx, ny=ny, nz=nz, origin=origin, A=A, B=B, C=C, values=values)

                i = j + 1
                continue
            end

            i += 1
        end
    end

    primvec === nothing && error("No PRIMVEC found")
    isempty(species) && error("No PRIMCOORD found")
    datagrid === nothing && error("No DATAGRID_3D found")

    structure = XSFStructure(primvec, species, positions)
    return structure, datagrid
end

function true_cell_vectors(datagrid)
    nx, ny, nz = datagrid.nx, datagrid.ny, datagrid.nz
    A = collect(Tuple(datagrid.A))
    B = collect(Tuple(datagrid.B))
    C = collect(Tuple(datagrid.C))
    At = A .* (nx / (nx - 1))
    Bt = B .* (ny / (ny - 1))
    Ct = C .* (nz / (nz - 1))
    return At, Bt, Ct
end

function grid_point(datagrid, i::Int, j::Int, k::Int)
    # Wannier90 half-open convention: i,j,k are 1-based -> (i-1)/N
    nx, ny, nz = datagrid.nx, datagrid.ny, datagrid.nz
    u = (i-1) / nx
    v = (j-1) / ny
    w = (k-1) / nz
    o = Tuple(datagrid.origin)
    At, Bt, Ct = true_cell_vectors(datagrid)
    return (
        o[1] + u*At[1] + v*Bt[1] + w*Ct[1],
        o[2] + u*At[2] + v*Bt[2] + w*Ct[2],
        o[3] + u*At[3] + v*Bt[3] + w*Ct[3],
    )
end

function _pad_periodic_2d(vals::AbstractMatrix)
    nx, ny = size(vals)
    out = Matrix{eltype(vals)}(undef, nx + 1, ny + 1)
    out[1:nx, 1:ny] = vals
    out[end, 1:ny] = vals[1, :]
    out[1:nx, end] = vals[:, 1]
    out[end, end] = vals[1, 1]
    return out
end

function _bilinear_periodic_2d(vals::AbstractMatrix{<:Real}, u::Real, v::Real)
    # vals is (nx, ny) sampled on u=(i-1)/nx, v=(j-1)/ny, periodic in both directions.
    nx, ny = size(vals)
    uu = u - floor(u)
    vv = v - floor(v)

    x = uu * nx + 1
    y = vv * ny + 1
    i0 = floor(Int, x)
    j0 = floor(Int, y)
    fx = x - i0
    fy = y - j0

    i1 = i0 + 1
    j1 = j0 + 1
    i0 = mod1(i0, nx)
    i1 = mod1(i1, nx)
    j0 = mod1(j0, ny)
    j1 = mod1(j1, ny)

    v00 = float(vals[i0, j0])
    v10 = float(vals[i1, j0])
    v01 = float(vals[i0, j1])
    v11 = float(vals[i1, j1])
    return (1 - fx) * (1 - fy) * v00 + fx * (1 - fy) * v10 + (1 - fx) * fy * v01 + fx * fy * v11
end

function slice_surface_data(dg, k::Integer; pad_periodic::Bool=true)
    # Returns X, Y, Z, V suitable for Makie.surface (works for non-orthogonal A,B).
    nx, ny, nz = dg.nx, dg.ny, dg.nz
    @assert 1 <= k <= nz

    O = StaticArrays.SVector{3,Float64}(Tuple(dg.origin))
    At, Bt, Ct = true_cell_vectors(dg)
    Atv = StaticArrays.SVector{3,Float64}(At...)
    Btv = StaticArrays.SVector{3,Float64}(Bt...)
    Ctv = StaticArrays.SVector{3,Float64}(Ct...)
    w = (k - 1) / nz
    p0 = O + w * Ctv

    X = Matrix{Float64}(undef, nx, ny)
    Y = Matrix{Float64}(undef, nx, ny)
    Z = Matrix{Float64}(undef, nx, ny)
    for j in 1:ny
        v = (j - 1) / ny
        for i in 1:nx
            u = (i - 1) / nx
            p = p0 + u * Atv + v * Btv
            X[i, j] = p[1]
            Y[i, j] = p[2]
            Z[i, j] = p[3]
        end
    end

    V = dg.values[:, :, k]
    if pad_periodic
        return _pad_periodic_2d(X), _pad_periodic_2d(Y), _pad_periodic_2d(Z), _pad_periodic_2d(V)
    end
    return X, Y, Z, V
end

function slice_heatmap_data(dg, k::Integer; ns::Integer=300, nt::Integer=300)
    # Resamples the (u,v) slice at index k onto an orthonormal in-plane (s,t) grid.
    # Output is suitable for Makie.heatmap(s, t, Vst).
    nx, ny, nz = dg.nx, dg.ny, dg.nz
    @assert 1 <= k <= nz

    At, Bt, _Ct = true_cell_vectors(dg)
    Atv = StaticArrays.SVector{3,Float64}(At...)
    Btv = StaticArrays.SVector{3,Float64}(Bt...)

    e1 = Atv / norm(Atv)
    Bt_perp = Btv - (dot(Btv, e1)) * e1
    nbp = norm(Bt_perp)
    nbp > 0 || error("A and B are colinear; cannot form an in-plane orthonormal basis")
    e2 = Bt_perp / nbp

    a11 = dot(Atv, e1)                 # = |A|
    a12 = dot(Btv, e1)                 # shear
    a22 = dot(Btv, e2)                 # = |B_perp|

    smin = min(0.0, a12)
    smax = a11 + max(0.0, a12)
    tmin = 0.0
    tmax = a22

    s = range(smin, smax; length=ns)
    t = range(tmin, tmax; length=nt)

    Vuv = dg.values[:, :, k]
    Vst = Matrix{Float64}(undef, ns, nt)

    for (it, tt) in enumerate(t)
        v = tt / a22
        for (is, ss) in enumerate(s)
            u = (ss - v * a12) / a11
            if 0.0 <= u < 1.0 && 0.0 <= v < 1.0
                Vst[is, it] = _bilinear_periodic_2d(Vuv, u, v)
            else
                Vst[is, it] = NaN
            end
        end
    end

    return collect(s), collect(t), Vst
end