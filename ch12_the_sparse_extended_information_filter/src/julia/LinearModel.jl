module LinearModel

export jacobian_motion, equation_motion,
jacobian_measurement, equation_measurement, inverse_measurement

H = 0.74
L = 2.83
a = 0.95 + L
b = 0.50

function jacobian_motion(θ,v,α,dt)
    [[1 0 -dt*(v*sin(θ)+v/L*tan(α)*(a*cos(θ)-b*sin(θ)))];
    [0 1 dt*(v*cos(θ)-v/L*tan(α)*(a*sin(θ)+b*cos(θ)))];
    [0 0 1]]
end

function equation_motion(θ,v,α,dt)
    v /= 1-tan(α)*H/L
    [dt*(v*cos(θ)-v/L*tan(α)*(a*sin(θ)+b*cos(θ))),
    dt*(v*sin(θ)+v/L*tan(α)*(a*cos(θ)-b*sin(θ))),
    dt*v/L*tan(α)]
end

function jacobian_measurement(μ, ind)
    n = length(ind)

    x = μ[2*ind .+ 2]
    y = μ[2*ind .+ 3]

    δ = [x y] .- μ[1:2]'
    q = sum(δ.^2,dims=2)

    J = [1 0;0 -1;0 1;1 0]*δ'
    J = [-J;zeros(1,n);-ones(1,n);J]
    J[1:2:end,:] .*= 1 ./ sqrt.(q')
    J[[2,4,8,10],:] .*= 1 ./ q'

    return reshape(J,(2,5,n))
end

function equation_measurement(μ,ind)
    n = length(ind)

    x = μ[2*ind .+ 2]
    y = μ[2*ind .+ 3]

    δ = [x y] .- μ[1:2]'
    q = sqrt.(sum(δ.^2,dims=2))

    [q atan.(δ[:,2],δ[:,1]).-μ[3].+π/2]'
end

function inverse_measurement(μ,z)
    c = cos.(μ[3] .+ z[2,:])
    s = sin.(μ[3] .+ z[2,:])

    return μ[1:2] .+ [(z[1,:].*s)'; (-z[1,:].*c)']
end

end
