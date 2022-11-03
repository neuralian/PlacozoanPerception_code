set_theme!(
    rowgap = 5,
    colgap = 5,
    fontsize = 16, 
    backgroundcolor = "#212946",
    Axis = (
        titlecolor = :white,
        xgridcolor = :white, 
        ygridcolor = :white, 
        backgroundcolor = "#212946",
        xgridstyle=:dash, 
        ygridstyle=:dash, 
        xgridwidth=0.35, 
        ygridwidth=0.35, 
        xtickalign=0, 
        ytickalign=0, 
        xtickcolor = :white, 
        ytickcolor = :white,
        bottomspinecolor = :white, 
        topspinecolor =:white,
        leftspinecolor =:white,
        rightspinecolor =:white,
        xlabelcolor=:white, 
        ylabelcolor = :white,
        xticklabelcolor = :white, 
        yticklabelcolor = :white,
        xgridvisible = true, 
        ygridvisible = true),
    Colorbar = (
        topspinevisible = false, 
        rightspinevisible = false,
        bottomspinevisible = false, 
        leftspinevisible = false,
        width = 12, 
        height = Relative(4/4), 
        labelcolor = :white,
        ticklabelcolor = :white,
        tickalign = 1, 
        tickcolor = :white, 
        labelpadding = -5, 
        ticklabelsize = 16),
    Legend = (
        bgcolor = ("#748AA6FF", 0.25), 
        labelcolor = :white, 
        titlecolor = :white, 
        framecolor=:white,
        tellheight = false, 
        tellwidth = false, 
        halign = :left, 
        valign = :bottom, 
        labelsize = 14, 
        linewidth = 2,
        margin = (10, 10, 10, 10),
        framevisible = true      
        )
)