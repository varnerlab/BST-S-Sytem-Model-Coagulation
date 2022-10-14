# include -
include("Include.jl")

# define a rect -
rect(w, h, x, y) = Shape(x .+ [0, w, w, 0, 0], y .+ [0, 0, h, h, 0])
gray(pct) = RGB(pct, pct, pct)

# generate some test data -
D = rand(10,10)

# what is the size of the data set -
(NR,NC) = size(D)
const marker_height = 1.0
const marker_width = 1.0

for j ∈ 1:NC

    y = (j-1)*(marker_height + 0.1)
    for i ∈ 1:NR
        
        x = (i-1)*(marker_width + 0.1)
        @show (x,y)

        shape = rect(marker_width, marker_height, x, y)
        plot!(shape, c=gray(i/255), label="")

    end

end
gui()

