using GLMakie
include("reader.jl")

file = "graphene_00002.xsf"

structure, dg = read_xsf(file)

fig = render_volume(dg)

k = 10

let
    s, t, Vst = slice_heatmap_data(dg, k; ns=300, nt=300)

    fig = Figure(size=(1000, 750))

    z = dg.origin[3] + k * dg.C[3] / 200

    ax = Axis(fig[1, 1], title="wannier function, k = $(k)", xlabel="s (Å)", ylabel="t (Å)", aspect=DataAspect())
    finite = filter(isfinite, vec(Vst))
    # max_val =  maximum(abs.(dg.values))
    max_val = maximum(abs.(finite))
    hm = heatmap!(ax, s, t, Vst; colormap=:balance, colorrange=(-max_val, max_val))

    # hm = heatmap!(ax, s, t, log10.(abs.(Vst) .+ 1e-16); colormap=:viridis)

    Colorbar(fig[1, 2], hm, label="value")

    save("graphene_00002_slice_k$(k).png", fig)
end

file = "graphene_00001.xsf"

structure, dg = read_xsf(file)

fig = render_volume(dg)

k = 10

let
    s, t, Vst = slice_heatmap_data(dg, k; ns=300, nt=300)

    fig = Figure(size=(1000, 750))

    z = dg.origin[3] + k * dg.C[3] / 200

    ax = Axis(fig[1, 1], title="wannier function, k = $(k)", xlabel="s (Å)", ylabel="t (Å)", aspect=DataAspect())
    finite = filter(isfinite, vec(Vst))
    # max_val =  maximum(abs.(dg.values))
    max_val = maximum(abs.(finite))
    hm = heatmap!(ax, s, t, Vst; colormap=:balance, colorrange=(-max_val, max_val))

    # hm = heatmap!(ax, s, t, log10.(abs.(Vst) .+ 1e-16); colormap=:viridis)

    Colorbar(fig[1, 2], hm, label="value")

    save("graphene_00001_slice_k$(k).png", fig)
end