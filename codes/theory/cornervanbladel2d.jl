# Predict corner singularity powers for 2D dielectric wedge
# using idea as in van Bladel's EM book, Sec 4.13.
# Barnett 11/6/25
#
# The "interior" material with relative permittivity epsilon
# occupies angle alpha. Exterior epsilon=1. Matching is
#     epsilon phi_n^- = phi_n^+
# where + denotes exterior and - interior side.
#
# For our SLP representation, the density power is the jump in normal
# derivative.

# Sorry no Project.toml for this. You'll need:
using Roots
using CairoMakie

al = pi/2   # corner angle
lam0 = 0.9    # a material ratio parameter (XG calls gamma I think)
epsfunc(lam) = (1+lam)/(1-lam)  # converts lanbda to epsilon perm interior
e = epsfunc(lam0)
@assert (e-1)/(e+1) ≈ lam0   # test inversion formula for eps from lam

"""
mismatch in shooting method for theta ODE for wedge of angle al,
relative premittivity e, with trial power g (so u''+g^2u = 0 in theta).
Even parity about the center of wedge for now, shoots to center of vacuum.
"""
function theta_shooting_even(al,e,g)
	# form is u(t) = cos(gt)                              for t < al/2,
	#                A cos(g(t-al/2)) + B sin(g(t-al/2))  for al/2 < t < pi.
	A = cos(g*al/2); B = -e*g*sin(g*al/2)     # val, deriv on right side of al/2
	if isinf(e) A=0.0; B=-1.0; end
	halfvac = pi-al/2                    # half the angle of vacuum exterior
	return -A*g*sin(g*halfvac) + B*g*cos(g*halfvac)  # deriv at t=pi (should be 0)
	# Note: one could equally return tan expression LHS-RHS of (4.66) van Bladel.
end
# *** todo: add theta_shooting_odd with other parity

nal = 200; als = collect(range(0,2pi,length=nal))    # how many angles
es = [1.0,2,5,10,38,50,100,Inf]    # list from van Bladel Fig 4.28.
gg = zeros(nal,length(es))
for (i,al) in enumerate(als)
	glast = 1.0
	for (j,e) in enumerate(es)
		#println("al=",al," e=",e)
		gg[i,j] = fzero(g -> theta_shooting_even(al, e, g), glast)
		glast = gg[i,j]      # use prev as init guess for power
		if gg[i,j]==0.0 gg[i,j]=NaN; end    # ignore if fzero locks on 0
	end
end
fig = Figure()
ax = Axis(fig[1,1],xlabel=L"wedge $\alpha$",ylabel=L"power $\gamma$",
		  title="Diel. wedge leading power (Fig. 4.28 from van Bladel EM book)")
ax.xticks = ([0,pi/2,pi,3pi/2,2pi], ["0","π/2","π","3π/2","2π"])
ax.limits = (nothing,(0,2))
for (j,e) in enumerate(es)
	lines!(als,gg[:,j], label=L"$\epsilon=%$e$")
end
axislegend(position=:rb)
display(fig)
save("fig4.28.svg",fig)

fig = Figure()    # reprod Xuanzhao's expt power plot
al = pi/2       # corner ang
nl = 100; lams = range(0.0,1.0,length=nl)
ne = 1    # nonzero powers (eigenvalues) to find in even or odd class; 1 for now
ax = Axis(fig[1,1],xlabel=L"ratio param. $\lambda$",ylabel=L"power $\gamma$",
		title=L"pot power-1 (=density power) for angle $\alpha=%$al$")
gg = zeros(nl,ne)
for (i,lam) in enumerate(lams)
	e = epsfunc(lam)
	println(lam,e)
	gg[i,1] = fzero(g -> theta_shooting_even(al, e, g), 1.0)
end
lines!(lams,gg[:,1] .- 1.0, color=:black)
display(fig)
save("powervslam.svg",fig)
