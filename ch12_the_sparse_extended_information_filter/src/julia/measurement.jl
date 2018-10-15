function measurement(z,c,μ,ξ,Ω,graph)
    n = size(graph,1)-1
    qnoise = diagm(0 => [5.0,0.02])
    nnew = maximum(c)-n

    if nnew > 0
        osize = length(ξ)
        nsize = osize+2*nnew
        Ω₁ = zeros(nsize,nsize)
        Ω₁[1:osize,1:osize] = Ω
        Ω = Ω₁
        append!(ξ,zeros(2*nnew))
        append!(μ,inverse_measurement(μ[1:3],z[:,c .> n])[:])
        padd = spzeros(Bool,nnew,n+1)
        padd1 = spzeros(Bool,n+nnew+1,nnew)
        graph = [[graph;padd] padd1]
    end

    jcb = jacobian_measurement(μ,c)
    ẑ = equation_measurement(μ,c)
    qnoiseinv = inv(qnoise)

    for (i,j) ∈ enumerate(c)
        h = jcb[:,:,i]
        dz = z[:,i]-ẑ[:,i]
        dz[2] = measure.(dz[2])
        ind = [1:3;2j+2;2j+3]

        ξ[ind] .+= h'*qnoiseinv*(dz+h*μ[ind])
        Ω[ind,ind] .+= h'*qnoiseinv*h
    end

    c = unique(c)
    rows = ones(length(c))
    cols = c .+ 1
    sg = size(graph,1)
    adjacency = sparse(rows,cols,true,sg,sg)
    graph .|= adjacency
    graph .|= adjacency'

    return Ω,graph

end
