include("Include.jl")

# what training index?
index = 10

# go -
(p, Yₘ, Y) = learn_optim(index; pₒ = nothing)