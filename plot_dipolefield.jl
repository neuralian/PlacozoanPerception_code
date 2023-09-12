# plot dipole field

fig = Figure(resolution = (800, 600))
ax = Axis(fig[1,1])
#lines!(predator.field)


# r = 10.:.01:100.
# E = [dipoleFieldstrength(i) for i in r]
# # lines!(r, E)

#  vv = cumsum(E[end:-1:1])[end:-1:1]

# lines!(r,vv)

# pov = pOpen(vv)
# lines!(r, pov*maximum(vv))

r = 1:400

pvolts = cumsum(predator.field[end:-1:1])[end:-1:1]
lines!(r, pvolts)

pof = pOpen(pvolts)
lines!(r, pof*maximum(pvolts))

E = [dipoleFieldstrength(i) for i in float(r)]
Evolts = cumsum(E[end:-1:1])[end:-1:1]
Evolts = Evolts*maximum(pvolts)/maximum(Evolts)  # scale to plot
lines!(r, Evolts)

display(fig)