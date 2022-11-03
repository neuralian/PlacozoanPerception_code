# Plot field and open state probability vs distance

include("makieTheme1.jl") 

#fontsize_theme = Theme(fontsize = 20)
set_theme!(Theme(fontsize = 24))

figF1 =Figure(resolution = (800,600))
axF1 = Axis(figF1[1,1], 
     yticks = vec(collect(0.0:250.0:750.0)), 
    title = "Field Strength and Likelihood of Range",
    xlabel = "Range (μm)",
    ylabel = "μv/cm")
    ylims!(-25.0, 1050.0)
axF1r = Axis(figF1[1,1], yaxisposition = :right,
   yticks = vec(collect(0.0:0.25:1.0)), 
   xlabel = "Range (μm)",
   ylabel = "Open State Probability")
ylims!(-0.025, 1.05)



    lines!(axF1, predator.field[1:175], color = "#2d82ba", linewidth=2.5 )
    lines!(axF1r, pOpenGivenFieldstrength(predator.potential[1:175].*1.0e-6),
    color =  "#a34003", linewidth=2.5  )
    lines!(axF1r, 1.0.-pOpenGivenFieldstrength(predator.potential[1:175].*1.0e-6),
    color =  "#c4a233" , linewidth=2.5  )

 


    # figF1[1,1] = axF1
     display(figF1)
     save("Field_fig.png", figF1, px_per_unit = 3 )
         #xticks = dticks,