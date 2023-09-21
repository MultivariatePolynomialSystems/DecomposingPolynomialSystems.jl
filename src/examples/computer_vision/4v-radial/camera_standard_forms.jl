include("../../../src/utils/utils.jl")

c = randn(ComplexF64, 3)
t = randn(ComplexF64, 2)
P = ct2sPrad(c,t)

P[1,1:3]*P[2,1:3]

# Uncalibrated form, Equation (5)
Ps = [randn(ComplexF64, 2, 4) for i=1:4]
A = vcat([transpose(Ps[i][1,:]) for i=1:4]...)
d = zeros(ComplexF64, 4)
for i = 1:4
    As = copy(A)
    As[i,:] = (1/Ps[1][1,1])*Ps[1][2,:]
    d[i] = det(As)
end
d, diagm(d)
# d = d/norm(d)
H = inv(diagm(d)*A)
H*diagm(d)*A
# det(H)
Ps = [Ps[i]*H for i=1:4]
Ps[1]


# Calibrated form, Equation (10)
cs = [randn(ComplexF64, 3) for i=1:4]
ts = [randn(ComplexF64, 2) for i=1:4]
Ps = [ct2Prad(cs[i], ts[i]) for i=1:4]
r1, r2 = Ps[1][1,1:3], Ps[1][2,1:3]
R1 = [r1 r2 cross(r1,r2)]
t1 = [-ts[1]; 0]
H1 = [R1 t1; 0 0 0 1]
P2H1 = Ps[2]*H1
t2 = -P2H1[1,4]*P2H1[1,1:3]
H2 = [P2H1[2,4]*eye(3) t2; 0 0 0 1]
H = H1*H2

[Ps[i]*H for i=1:4]

# Upright case
