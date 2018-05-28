%simulation of random variable (x+cos?, y+sin?)
function r = simmov(s1, s2, n)
    theta = normrnd(0, sqrt(s2), [n,1])
    xy = mvnrnd([0 0], s1*eye(2), n)
    r = xy + [cos(theta) sin(theta)]
end
