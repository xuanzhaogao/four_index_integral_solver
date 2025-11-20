using LinearAlgebra, IterativeSolvers
using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie
using CSV, DataFrames

rects = [BI.square(-1.0, -0.5), BI.square(0.0, -0.5)]

n_panel = 8
n_quad = 16

n_adapts = 10:5:20
epses = range(1.1, 300.0, length = 100)

df = CSV.write("data/condition_number_threaded.csv", DataFrame(n_adapt = [], gamma = [], eps = [], cond = [], niter = [], res = []))

for n_adapt in n_adapts

    gammas = zeros(length(epses))
    conds = zeros(length(epses))
    niters = zeros(length(epses)) 
    reses = zeros(length(epses))

    Threads.@threads for i in eachindex(epses)
        eps_box = epses[i]
        gamma = (eps_box - 1.0) / (eps_box + 1.0)
        gammas[i] = gamma

        eps_boxes = [2.0, eps_box]
        eps_src = 2.0
        r_src = (-0.1, 0.4)

        mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panel, n_quad, n_adapt)
        lhs = BI.Lhs_dielectric_mbox2d(mbox)
        rhs = BI.Rhs_dielectric_mbox2d(mbox, r_src, eps_src)

        cond_number = cond(lhs)

        # try gmres to solve the system
        sigma, history = gmres(lhs, rhs, reltol = 1e-10, log = true, restart = 100)

        println("n_adapt = $(n_adapt), gamma = $(gamma), eps = $(eps_box), cond = $(cond_number), niter = $(length(history.data[:resnorm])), res = $(history.data[:resnorm][end])")

        conds[i] = cond_number
        niters[i] = length(history.data[:resnorm])
        reses[i] = history.data[:resnorm][end]
    end

    CSV.write("data/condition_number_threaded.csv", DataFrame(n_adapt = n_adapt, gamma = gammas, eps = epses, cond = conds, niter = niters, res = reses), append = true)

end