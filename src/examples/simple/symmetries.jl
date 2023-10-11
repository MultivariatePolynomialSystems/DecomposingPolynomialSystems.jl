# using DecomposingPolynomialSystems

# #-----------------------------------------------------------------#

# @var x a b
# F = System([x^2 + a*x + b]; parameters = [a,b])

# symmetries, Fs = symmetries_fixing_parameters(F; degree=1, param_dep=false)

# #-----------------------------------------------------------------#

# @var x y p
# F = System([x^2 + x + p, x + y + p]; parameters = [p])

# symmetries, Fs = compute_symmetries(F, degree=1, std_form=true)
# symmetries.basis
# symmetries.coefficients

# symmetries = compute_symmetries(F, degree=1)
# symmetries

# #-----------------------------------------------------------------#

# @var x a b
# F = System([(x-a)^4+a*(x-a)^3+b*(x-a)^2+a*(x-a)+1]; parameters=[a,b])

#-----------------------------------------------------------------#

using DecomposingPolynomialSystems
@var x y z w a b c
eqs = [2*w^2+1, 2*x+4*w*y+2*a, -3*y^2-z^2+4*a*w*y+2*b, -w*y^3-3*w*y*z^2-a*y^2-a*z^2+2*b*w*y+2*c]
F = System(eqs; variables=[x,y,z,w], parameters=[a,b,c])

# some scalings jump from one component to another => problem in find_solution_id
scalings = scaling_symmetries(F)
scalings.grading[2][2]

F = run_monodromy(F)

symmetries_fixing_parameters!(F; degree_bound=1, param_dep=false)

gr = scaling_symmetries(F)
gr.U[2]
gr.U[1]

A = [0 -0.707106781186548*im -0.707106781186548*im; 
     0.707106781186547*im 0.5 -0.5;
     -0.707106781186548*im 0.5 -0.5]
B = [-1 0 0; 0 -1 0; 0 0 1]
A*B-B*A


# [w, x, y, z]
# [w, -sqrt(-0.5)*(y+z), -0.5*sqrt(-0.5)*(a-x) - 0.5*z, sqrt(-0.5)*a]

for i in eachindex(symmetries)
    println(symmetries[i])
    println()
end


symmetries, _ = compute_symmetries(F, xp0, degree=3)
symmetries


