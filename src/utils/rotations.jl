export q2sR, q2R, c2sR, c2R, q2c, R2c

function q2sR(q)
    a,b,c,d = q
    return [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
            2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
            2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2]
end

function q2R(q)
    a,b,c,d = q
    return 1/(a^2+b^2+c^2+d^2)*q2sR(q)
end

function c2sR(c)
    x, y, z = c
    return [1+x^2-y^2-z^2 2*(x*y-z) 2*(y+x*z);
            2*(x*y+z) 1-x^2+y^2-z^2 2*(y*z-x);
            2*(x*z-y) 2*(x+y*z) 1-x^2-y^2+z^2]
end

function c2R(c)
    x, y, z = c
    return 1/(1+x^2+y^2+z^2)*c2sR(c)
end

q2c(q) = [q[2]/q[1], q[3]/q[1], q[4]/q[1]]
R2c(R) = 1/(1+tr(R))*xx2v((R-transpose(R)))