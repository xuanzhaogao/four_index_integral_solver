using BoundaryIntegral
import BoundaryIntegral as BI
using JLD2
using CairoMakie

res = load(joinpath(@__DIR__, "data/single_box3d_ref.jld"))["res_ref"]
res_2 = load(joinpath(@__DIR__, "data/single_box3d_4_12_6_8.jld"))["res"]
trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

X = trgs[1, :]
Y = trgs[2, :]
Z = trgs[3, :]

F = res[4]

begin
    θ = range(-π/2, π/2; length=100)
    φ = range(-π, π; length=100)

    Θ = [t for t in θ, p in φ]
    Φ = [p for t in θ, p in φ]

    # 3. 判断哪些经纬度点属于“前半球”（从 +x 方向看）
    #    球面坐标 -> 笛卡尔: (x,y,z) = (cosθ cosφ, cosθ sinφ, sinθ)
    #    从 +x 看：前半球当且仅当 x >= 0
    mask_front = @. cos(Θ) * cos(Φ) >= 0

    # 4. 构造一个带 alpha 通道的颜色矩阵
    #    可见部分 alpha=1，不可见部分 alpha=0
    colors = RGBAf.(F .- minimum(F), 0, 0, 1)  # 这里举例只把数值编码放在 red，你可以换 colormap
    # 若要用 colormap，后面可直接用 heatmap 的 colormap 参数，不用手动 RGBA 映射。

    # 改为只改变透明度
    alpha_mat = map(x -> x ? 1.0f0 : 0.0f0, mask_front)
    # RGBA with alpha
    C = RGBAf.(F, F, F, alpha_mat)  # 实际颜色用 heatmap colormap 覆盖，此处 alpha 保存即可

    # 5. 画图：先画 heatmap，再盖一个圆形轮廓
    fig = Figure(resolution=(700,700))
    ax = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")

    # 把 (φ, θ) 映射到一个单位圆上展示
    # 最简单的做法：直接用 (x = cosθ cosφ, y = cosθ sinφ)
    X = @. cos(Θ) * cos(Φ)
    Y = @. cos(Θ) * sin(Φ)

    heatmap!(ax, X, Y, F;
            colormap=:viridis,
            transparency=true,
            alpha=alpha_mat)

    # 边界圆
    t = range(0, 2π; length=400)
    lines!(ax, cos.(t), sin.(t), color=:black)

    xlims!(ax, -1.05, 1.05)
    ylims!(ax, -1.05, 1.05)
    hidedecorations!(ax)
    hidespines!(ax)
end