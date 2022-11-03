
using GLMakie

N = 2 # default number of points in circle decomposition
M = 2   # number of cells to draw, must == length(Φ)
A = fill(fill(Point2f0(0,0), N), M)
B    = fill(fill(Point2f0(0,0), N), M)

for i in 1:M

   # each element of A is a vector of 64 points
   A[i] = decompose(Point2f0, Circle( Point2f0(rand(),rand()), 1.0) )[1:N]

   for j in 1:N
      # copy points 1 by 1
      B[i][j]= A[i][j]
   end

end
