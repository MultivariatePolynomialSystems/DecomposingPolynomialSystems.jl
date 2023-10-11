using DecomposingPolynomialSystems, JLD2

@var x[1:3] y[1:3] z[1:3] t[1:2,1:2] α[1:13,1:4] X[1:13] l[1:2,1:13,1:4]

cs = vcat([[0,0,0]], [[x[i], y[i], z[i]] for i=1:3])
ts = [[0;0], [1;0], t[:,1], t[:,2]]
Ps = [ct2sPrad(cs[i], ts[i]) for i = 1:4]
eqs = Vector{Expression}([])
for i in 1:13
    for j in 2:4
        append!(eqs, α[i,j]*l[:,i,j] - Ps[j]*[α[i,1]*l[:,i,1]; X[i]; 1])
    end
end

F = load_object("src/examples/vision/4v-radial/4v-radial.jld2")
F.system = System(eqs; variables=vcat(M2V([x y z]'), M2V(t), M2V(α), M2V(X)), parameters=M2V(l))

deck = symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false, logging=true)

show_map(deck, 12)

for a in 1:3584
    b = F.deck_permutations[12][a]
    sol1 = F.samples.solutions[:,a,1]
    sol2 = F.samples.solutions[:,b,1]
    x₁ = sol1[findfirst(x->x==α[6,3], unknowns(F))]
    x₂ = sol2[findfirst(x->x==α[6,3], unknowns(F))]
    if abs(x₁-x₂) > 1e-5
        println("FOUND")
    end
end