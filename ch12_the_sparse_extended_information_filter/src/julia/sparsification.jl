function sparsification(m₀,m₁,μ,ξ,Ω,graph)

    sg = size(graph,1)
    rows = ones(length(m₁))
    cols = m₁ .+ 1
    adjacency = sparse(rows,cols,true,sg,sg)
    graph .⊻= adjacency .| adjacency'

    edges = [[i+1 j+1] for i ∈ union(m₀,m₁),j ∈ m₁ if i ≠ j]
    edges = vcat(edges...)
    adjacency = sparse(edges[:,1],edges[:,2],true,sg,sg)
    graph .|= adjacency .| adjacency'

    ind₀ = [[2i+2,2i+3] for i ∈ sort(m₀)]
    ind₀ = vcat(ind₀...)
    ind₁ = [[2i+2,2i+3] for i ∈ sort(m₁)]
    ind₁ = vcat(ind₁...)

    Ω₂ = Ω[[1:3;ind₀;ind₁],ind₁]
    l₁ = Ω₂ / Symmetric(Ω[ind₁,ind₁]) * Ω₂'

    Ω₃ = Ω[[1:3;ind₀;ind₁],[1:3;ind₁]]
    l₂ = Ω₃ / Symmetric(Ω[[1:3;ind₁],[1:3;ind₁]]) * Ω₃'

    Ω₄ = Ω[:,1:3]
    l₃ = Ω₄ / Symmetric(Ω[1:3,1:3]) * Ω₄'

    Ω₁ = copy(Ω)
    Ω[[1:3;ind₀;ind₁],[1:3;ind₀;ind₁]] .-= l₁
    Ω[[1:3;ind₀;ind₁],[1:3;ind₀;ind₁]] .+= l₂
    Ω .-= l₃

    ξ .+= (Ω-Ω₁)*μ
end
