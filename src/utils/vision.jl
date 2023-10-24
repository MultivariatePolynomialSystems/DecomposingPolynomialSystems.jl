export ct2P, ct2sP, ct2Prad, ct2sPrad

ct2P(c,t) = [c2R(c) t]
ct2sP(c,t) = (1+c[1]^2+c[2]^2+c[3]^2)*ct2P(c,t)
ct2Prad(c,t) = [c2R(c)[1:2,:] t]
ct2sPrad(c,t) = (1+c[1]^2+c[2]^2+c[3]^2)*ct2Prad(c,t)
