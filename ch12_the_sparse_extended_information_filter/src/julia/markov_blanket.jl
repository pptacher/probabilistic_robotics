function markov_blanket(n,m₀,graph)
    sg = size(graph,1)
    mb1::Array{UInt16,1} = m₀
    mb2::Array{UInt16,1} = findall(graph[n,:])

    mb12 = intersect(mb1,mb2,sg)
    if length(mb12) == 0
        mb = bfs(n,graph)
        if length(mb) > 0
            mb2 = union(mb2,mb,sg)
        end
    end

    mb2 = [1;n;mb2]
    union(mb1,mb2,sg)
end

function Base.intersect(a::AbstractArray{UInt16,1},b::AbstractArray{UInt16,1},s::Int)
    x = zeros(Int,s)
    y = zeros(Int,s)
    x[a] .= 1
    y[b] .= 1
    findall((x .* y) .== 1)
end

function Base.union(a::AbstractArray{UInt16,1},b::AbstractArray{UInt16,1},s::Int)
    x = zeros(Int,s)
    y = zeros(Int,s)
    x[a] .= 1
    y[b] .= 1
    findall((x .+ y) .> 0)
end

function bfs(n,graph)
    sg = size(graph,1)
    d = -ones(Int,sg)
    d[1] = 0
    q = zeros(Int,sg)
    i = 1
    j = 2
    q[i] = 1
    found = false
    p = zeros(Int,0)
    while  !found && i < j
        curr = q[i]
        i += 1
        next = findall(graph[n,:] .* (d.==-1))
        s1 = length(next)
        q[j:j+s1-1] = next
        j += s1
        d[next] .= curr
        found=d[n]==curr
    end

    if found
        p = curr
        prev = d[curr]
        while prev ≠ 1
            p = [p prev]
            prev = d[prev]
        end
    end
    return p
end
