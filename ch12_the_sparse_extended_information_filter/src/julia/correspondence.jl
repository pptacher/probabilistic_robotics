function correspondence(z,m₀,μ,Ω,graph)
    k = size(z,2)
    if k == 0
        return zeros(0)
    end
    n = (length(μ)-3)÷2
    qnoise = diagm( 0 => [5.0,0.02])
    co = zeros(Int,k)

    if n > 0
        t = μ[4:end]
        t₁ = t[1:2:end]
        t₂ = t[2:2:end]
        lst = findall((abs.(t₁.-μ[1]) .≤ 75) .*
                (abs.(t₂.-μ[2]) .≤ 75) .*
                ((t₁.-μ[1])*cos(μ[3]) .+ (t₂.-μ[2])*sin(μ[3]) .≥ -5))
        a = length(lst)
        jcb = jacobian_measurement(μ,lst)
        ẑ = reshape(equation_measurement(μ,lst),(2,1,a))

        q = SharedArray{Float64,2}((2a,2))

        @sync begin
            for p in procs(q)
                @async remotecall_wait(compute_cov_shared!,p,
                                        Ω,lst,m₀,graph,q,jcb,qnoise)
            end
        end

        f = reshape(permutedims(ẑ.-z,[1,3,2]),(2*a,k))
        f[2:2:end,:] = measure.(f[2:2:end,:])
        d₁ = zeros(2a,k)
        cc, ia = intersect_ind(lst,m₀,n)

        @inbounds for j ∈ ia
            d₁[2j-1:2j,:] = q[2j-1:2j,:]*f[2j-1:2j,:]
        end

        ind = [[2j-1,2j] for j ∈ ia]
        d₂ = d₁[vcat(ind...),:]
        d₂ .= d₂.^2
        d₃ = [d₂[i,j]+d₂[i+1,j] for i ∈ 1:2:2length(ia),j ∈ 1:k]

        co = argmin(d₃,dims=1)
        co = co[:]
        co = getindex.(co,1)
        m₁ = minimum(d₃,dims=1)[:]
        c₁ = cquantile(Chisq(2), .50)
        c₂ = cquantile(Chisq(2), .05)
        newl = findall(m₁ .> c₁)
        new₁ = findall(m₁[newl] .< c₂)
        ia = ia[setdiff(1:length(ia),co[newl[new₁]])]
        co = cc[co]

        if (length(newl) > 0 && length(cc) < a) || length(cc) == 0
                if length(cc) == 0
                        newl = 1:k
                end
                ib = setdiff(1:a,ia)
                lst = lst[ib]
                @inbounds for j ∈ ib
                        d₁[2j-1:2j,newl] = q[2j-1:2j,:]*f[2j-1:2j,newl]
                end
                if length(ib) > 0
                        ind = [[2j-1,2j] for j ∈ ib]
                        d₂ = d₁[vcat(ind...),newl]
                        d₂ = d₂.^2
                        d₃ = [d₂[i,j]+d₂[i+1,j] for i ∈ 1:2:2*length(ib),j ∈ 1:length(newl)]

                        co₁ = argmin(d₃,dims=1)
                        co₁ = co₁[:]
                        co[newl] = lst[getindex.(co₁,1)]
                        m₁ = minimum(d₃,dims=1)[:]
                        new₁ = []
                        newl = newl[m₁ .> c₂]
                end
        end
    else
        newl = 1:k
        new₁ = []
    end

    if length(new₁) > 0
            newl = newl[[i for i ∈ 1:length(newl) if i ∉ new₁]]
    end
    if length(newl) > 0
            co[newl] = n .+ [i for i ∈ 1:length(newl)]
    end
    return co
end

@everywhere include("markov_blanket.jl")

@everywhere function myrange(q)
        idx = indexpids(q)
        if idx == 0
            return 1:0
        end
        nchunks = length(procs(q))
        splits = [round(Int,s) for s ∈ range(0,stop=size(q,1)÷2,length=nchunks+1)]
        splits[idx]+1:splits[idx+1]
end

@everywhere function compute_cov!(Ω,lst,m₀,graph,q,jcb,qnoise,irange)
        @inbounds for i ∈ irange
            ind = markov_blanket(lst[i]+1,m₀.+1,graph)
            ind = ind[2:end] .- 1
            ind = [[2j+2,2j+3] for j ∈ ind ]
            ind = vcat(ind...)
            ind1 = findfirst(ind.==2*lst[i]+2)
            f = fxm(ind1,3+length(ind))
            j = jcb[:,:,i]
            #drastic improvements with symmetric matrices specific solver.
            s = j*f / Symmetric(Ω[[1:3;ind],[1:3;ind]]) *f'*j'
            vals,vec = eigen(Symmetric(s+qnoise))
            q[2i-1:2i,:] = diagm( 0 => sqrt.(1 ./ vals) )*vec'
        end
end

@everywhere function compute_cov_shared!(Ω,lst,m₀,graph,q,jcb,qnoise)
        compute_cov!(Ω,lst,m₀,graph,q,jcb,qnoise,myrange(q))
end

@everywhere function fxm(p,n)
        sparse(1:5,[1,2,3,3+p,4+p],1,5,n)
end

function measure(θ)
    tmp = θ%2π
    tmp-min(tmp÷π,1)*2π
end

function intersect_ind(a,b,n)
    #assume a is sorted vectors of positive integers
    u = zeros(Int,2,n)
    v = zeros(Int,n)
    u[1,a] .= 1
    u[2,a] = 1:length(a)
    v[b] .= 1
    ps = u[1,:] .* v
    res =  findall(ps .≠ 0)
    res, u[2,(ps .* u[2,:]) .> 0]
end
