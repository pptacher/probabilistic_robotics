function motion(v,α,dt,m₀,μ,ξ,Ω,graph)

    n = length(ξ)
    rnoise = diagm(0 => [8e-4,8e-4,1e-6])
    jcb = jacobian_motion(μ[3],v,α,dt)

    ψ = inv(jcb)-I
    ind = [[2*i+2;2*i+3] for i ∈ sort(m₀)]
    ind = vcat(ind...)
    ind = [1:3;ind]
    Ω₁ = Ω[ind,1:3]
    ψ₁ = Ω₁*ψ

    λ = zeros(n,n)
    λ[ind,1:3] = ψ₁
    λ[1:3,ind] .+= ψ₁'
    λ[1:3,1:3] .+= ψ'*Ω[1:3,1:3]*ψ
    Φ = Ω + λ

    Φₓ = Φ[1:3,ind]
    κ = zeros(n,n)
    κ[ind,ind] = Φₓ'*inv(inv(rnoise)+Φ[1:3,1:3])*Φₓ
    Ω .+= λ - κ

    δ = equation_motion(μ[3],v,α,dt)
    ξ[ind] .+= Ω[ind,1:3]*δ + (λ-κ)[ind,:]*μ
    μ[1:3] += δ

    if length(m₀)>0
        edges = [[i+1 j+1] for i ∈ m₀,j ∈ m₀ if i ≠ j]
        edges = vcat(edges...)
        if size(edges,1)>0
            sg = size(graph,1)
            adjacency = sparse(edges[:,1],edges[:,2],true,sg,sg)
            graph .|= adjacency
        end
    end
end
