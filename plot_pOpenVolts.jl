using GLMakie

fig = Figure(resolution = (800,600))
ax = Axis(fig[1,1])

eMax = 6.0e4
# xlims!(-2.0e-8, eMax)
# ylims!(0.0, 1.0)
de = 1.0e2
e = -(eMax/2.0):de:eMax 
lines!(e, pOpen(e))

lines!([0., eMax], 0.1*[1,1])
lines!(physics.σ*[1,1], [0.,1.])

display(fig)