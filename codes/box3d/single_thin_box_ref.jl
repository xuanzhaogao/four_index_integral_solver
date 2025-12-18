include(joinpath(@__DIR__, "single_box3d_utils.jl"))

# consider three sets of points for accuracy check:
# plan z = 2.0
n = 100

trgs_1 = zeros(3, n * n)
trgs_2 = zeros(3, n * n)
trgs_3 = zeros(3, n * n)

x = range(-1.0, 1.0, length = n)
for i in 1:n
    for j in 1:n
        trgs_1[:, (i - 1) * n + j] = [x[i], x[j], 2.0]
        trgs_2[:, (i - 1) * n + j] = [x[i], x[j], 0.6]
        trgs_3[:, (i - 1) * n + j] = [x[i], 0.0, x[j]]
    end
end


tbox, sigma, gi, n_val, n_iter = solve_single_thin_box3d(4.0, 5.0, 5.0, 1.0, 5, 5, 1, 4, 4, 10, 10, (0.2, 0.3, 0.4))

