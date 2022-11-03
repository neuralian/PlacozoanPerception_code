  
set_theme!(
  rowgap = 5,
  colgap = 5, 
  fontsize = 16, 
  font =:sans,
  Axis = (
      xgridcolor = :black, 
      ygridcolor = :black, 
      backgroundcolor = :white,
      xgridstyle=:dash, 
      ygridstyle=:dash, 
      xgridwidth=0.15, 
      ygridwidth=0.15, 
      xtickalign=1, 
      ytickalign=1, #xtickcolor = :red, ytickcolor = :white,
      xgridvisible = true, 
      ygridvisible = true),
  Colorbar = (
      topspinevisible = false,
      rightspinevisible = false,
      bottomspinevisible = false, 
      leftspinevisible = false,
      width = 12, 
      height = Relative(4/4), 
      tickalign = 1, 
      labelpadding = -5),
  Legend = (
      tellheight = false, 
      tellwidth = false, 
      halign = :left, 
      valign = :bottom, 
      labelsize = 14, 
      linewidth = 2,
      margin = (10, 10, 10, 10),
      framevisible = true)
)