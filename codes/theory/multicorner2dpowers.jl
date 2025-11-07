# Predict corner singularity power sequence for multijunction 2D dielectric
# generalized wedge, using idea as in van Bladel's EM book, Sec 4.13.
# Barnett 11/7/25, generalizing cornervanbladel2d.jl to use 2x2 det of ODE
# transmission matrix (or one element of such matrix when PEC).
#
# There are nm>=2 materials including vacuum (nm=2 is a plain diel corner).
# The last material is optionally PEC (eps=Inf).
# The relative permittivities are length-nm vector e, and angles a length nm-1
# vector (the last angle defined implicitly).
# Matching at each junction is
#     eps_j phi_n^- = eps_{j+1} phi_n^+
# where -(+) denotes  theta just below (above).
# For our SLP representation, the density power is the jump in normal
# derivative.

# Sorry no Project.toml for this. You'll need:
using Roots, LinearAlgebra
using CairoMakie

"""
	theta_ODE_det(a,e,g)

	return determinant for theta ODE transmission for multi-junction dielectric
	relative premittivity vector `e` and angles `a` (length one less than `e`,
	since the last angle is known), and trial power `g`
	(so u''+g^2u = 0 is the theta ODE in each material).
	Zero is returned iff `g` is a (nonlin) eigenvalue of the periodic ODE on
	[0,2pi].
	The last `e` value (only) may be `Inf`, in which case Dirichlet BCs are used
	on its angle endpoints (and the incoming u' to final u value map used).
	Note: Parity issues avoided due to solving the 2x2 transmission matrix
	on the entire [0,2pi] domain. (All parities are combined.)
"""
function theta_ODE_det(a::AbstractVector,e::AbstractVector,g)
	M = Matrix(1.0I,2,2)     # initial I transm matrix (at theta=0^+)
	aa = [a; 2pi-sum(a)]       # add the last angle
	@assert aa[end]>=0 "angles ($a) sum to bigger than 2pi!"
	nm = length(e)           # no. materials
	@assert length(aa)==nm "input vectors wrong lengths!"
	nexte = circshift(e,-1)   # list of j+1'th epsilons
	for (j,e) in enumerate(e)
		if !isinf(e)          # if last material conductor (u=0) skip it
			c = cos(g*aa[j]); s = sin(g*aa[j])      # where nonlin in g enters
			# update transm matrix by propagate by angle a then scale u'...
			soverg = g==0.0 ? aa[j] : s/g       # handle g=0 (u=at+b)
			M = [1.0 0; 0 e/nexte[j]] * [c soverg; -g*s c] * M
		end
	end
	return isinf(e[end]) ? M[1,2] : det(M-I)   # use u' -> u map if Dir BCs 
end

@assert abs( theta_ODE_det([3pi/2], [1.0,Inf], 2/3) )<1e-14  # PEC right-angle
# eps=2 right-angle...
@assert abs( theta_ODE_det([pi/2], [2.0,1.0], 1.1066007580762274) )<1e-14
g1 = fzero(g->theta_ODE_det([pi/2], [2.0,1.0], g),1.0)

gg = 0:0.01:3    # also used below
if false # plot a sweep over g to see higher zeros (Note tuples stop the broadcasting):
fig,ax,p = lines(gg, theta_ODE_det.(([pi/2],), ([2.0,1.0],), gg))
ax.xlabel=L"power $\gamma$"; ax.ylabel="det"
display(fig)
end

# higher powers w/o rootfinding...
angs = range(0,2pi,length=200)
eps = [10.0,1.0]     # diel wedge (ang is for diel)
#eps = [1.0, Inf]   # PEC corner (ang is for vacuum)
dd = [theta_ODE_det([a], eps, g) for a in angs, g in gg]
fig,ax,p = heatmap(angs,gg,-log.(abs.(dd).+1e-10))    # log highlights the zeros
p.colormap=:jet; p.colorrange=(-10,10)
ax.xticks = ([0,pi/2,pi,3pi/2,2pi], ["0","π/2","π","3π/2","2π"])
display(fig)
save("diel_higher_powers.svg",fig)

nal = 200; als = collect(range(0,2pi,length=nal))    # how many angles
es = [1.0,2,5,10,38,50,100] #,Inf]    # list from van Bladel Fig 4.28.
gg = zeros(nal,length(es))
for (i,al) in enumerate(als)
	glast = 1.0
	for (j,e) in enumerate(es)
		#println("al=",al," e=",e)
		gg[i,j] = fzero(g -> theta_ODE_det([al], [e,1.0], g), glast-0.1)  # hack to lock root
		glast = gg[i,j]      # use prev as init guess for power
		if gg[i,j]==0.0 gg[i,j]=NaN; end    # ignore if fzero locks on 0
	end
end
fig = Figure()
ax = Axis(fig[1,1],xlabel=L"wedge $\alpha$",ylabel=L"power $\gamma$",
		  title="Diel. wedge leading power, both parities (Fig. 4.28 from van Bladel)")
ax.xticks = ([0,pi/2,pi,3pi/2,2pi], ["0","π/2","π","3π/2","2π"])
ax.limits = (nothing,(0,2))
for (j,e) in enumerate(es)
	lines!(als,gg[:,j], label=L"$\epsilon=%$e$")
end
axislegend(position=:rt)
display(fig)
save("fig4.28_bothparities.svg",fig)


epsfunc(lam) = (1+lam)/(1-lam)  # converts lanbda to epsilon perm interior
lam0 = 0.9    # a material ratio parameter (XG calls gamma I think)
e = epsfunc(lam0)
@assert (e-1)/(e+1) ≈ lam0   # test inversion formula for eps from lam

begin fig = Figure(size=(300,1000))    # reprod Xuanzhao's expt power plot
al = pi/2       # corner ang
nl = 100; lams = range(0.0,1.0,length=nl)
ne = 5    # nonzero powers (eigenvalues) to find in even or odd class; 1 for now
ax = Axis(fig[1,1],xlabel=L"ratio param. $\lambda$",ylabel=L"power $\gamma$",
		title=L"pot powers-1 (=density powers) for angle $\alpha=%$al$")
gg = zeros(nl,ne)
for (i,lam) in enumerate(lams)
	e = epsfunc(lam)
	gguess = [0.8,1.2,1.8,2.9,3.2]    # kinda just above and below odd integers... hack
	for n=1:ne
		gg[i,n] = fzero(g -> theta_ODE_det([2pi-al], [1.0,e], g), gguess[n])  # allows PEC case
	end
end
for n=1:ne lines!(lams,gg[:,n] .- 1.0, color=:black) end   # show all powers
display(fig)
save("dielrightangle_powers_vs_lam.svg",fig)
# redo using log det to check:
gg = 0:0.01:3  # powers
fig = Figure(size=(300,1000))
ax = Axis(fig[1,1],xlabel=L"ratio param. $\lambda_3$",ylabel=L"power $\gamma$",
		title=L"$$pot powers (log det) for right-ang triple junc")
dd = [theta_ODE_det([2pi-al], [1.0,eps], g) for eps in epsfunc.(lams), g in gg]
p = heatmap!(lams,gg,-log.(abs.(dd).+1e-10))    # log highlights the zeros
p.colormap=:jet; p.colorrange=(-10,10)
display(fig)
end

# *** could do Boyd rootfinding on the analytic func of g to get lowest few (real) powers

# multijunction test with halfspace of vacuum meeting fixed eps_2 and variable eps_3:
fig = Figure(size=(300,1000))
al = [pi, pi/2]       # corner angs: materials 1=vacuum,2=fixed, 3=variable eps
e2 = 2.0   # fixed 
nl = 100; lams3 = range(0.0,1.0,length=nl)   # for eps3
ne = 5    # nonzero powers (eigenvalues) to find in even or odd class; 1 for now
ax = Axis(fig[1,1],xlabel=L"ratio param. $\lambda_3$",ylabel=L"power $\gamma$",
		title=L"pot powers-1 (=density powers) for right-ang triple junc")
gg = zeros(nl,ne)
for (i,lam3) in enumerate(lams3)
	e3 = epsfunc(lam3)
	gguess = [0.8,1.2,1.8,2.9,3.2]    # kinda just above and below odd integers... hack
	for n=1:ne
		gg[i,n] = fzero(g -> theta_ODE_det(al, [1.0,e2,e3], g), gguess[n])  # allows PEC case
	end
end
for n=1:ne lines!(lams,gg[:,n] .- 1.0, color=:black) end   # show all powers
display(fig)
save("dieltripjunc_powers_vs_lam.svg",fig)
# redo using log det to check:
gg = 0:0.01:3  # powers
fig = Figure(size=(300,1000))
ax = Axis(fig[1,1],xlabel=L"ratio param. $\lambda_3$",ylabel=L"power $\gamma$",
		title=L"$$pot powers (log det) for right-ang triple junc")
dd = [theta_ODE_det(al, [1.0,e2,eps3], g) for eps3 in epsfunc.(lams3), g in gg]
p = heatmap!(lams3,gg,-log.(abs.(dd).+1e-10))    # log highlights the zeros
p.colormap=:jet; p.colorrange=(-10,10)
display(fig)
save("dieltripjunc_det_vs_lam.svg",fig)
end

