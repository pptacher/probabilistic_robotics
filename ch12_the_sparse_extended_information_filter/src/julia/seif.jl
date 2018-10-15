using LinearAlgebra
using LinearModel
using SparseArrays, SharedArrays
using HDF5
using Profile
using Distributed
using Distributions
using Plots
using Printf

macro print_array(q)
    return :(Base.print_array(stdout,round.( $(esc(q)) ,sigdigits=4));println())
end

include("seif_fct.jl")

function seif()
    filename1 = "aa3_lsr2.h5"
    filename2 = "aa3_dr.h5"
    filename3 = "z.h5"
    include_dir = "./data/"

    dt = 25e-3
    graph = spzeros(Bool,1,1)
    qnoise = diagm(0 => [5, 0.02])

    speed, steering, timeit = h5open(include_dir * filename2, "r") do file
          (read(file, "/speed"), read(file, "/steering"), read(file, "/time"))
    end

    timelsr = h5open(include_dir * filename1, "r") do file
          read(file, "/timeLsr")
    end
    llsr = size(timelsr,1)

    z = h5open(include_dir * filename3, "r") do file
          read(file, "/zz")
    end

    μ = zeros(3)
    ξ = zeros(3)
    Ω = 10e4*diagm(0 => [1,1,1])
    m₀ = zeros(Int,0)
    n = 20

    stindex = 1
    j = findfirst(timelsr .>= timeit[stindex])[1]

    gr(reuse=true)
    plt = plot([0;.1],[0;.2];axis=false,
    leg=false,
    aspect_ratio=1,
    grid=false)
    xlims!((-200,200))
    ylims!((-50,300))

    poses = zeros(3,5000)
    #ξ_hist = zeros(3,500)
    t₁ = time()

    for i ∈ stindex:size(timeit,1)
        t₂ = time()

        if i > stindex && i%5000==0
            plot!(poses[1,:],poses[2,:],line = (:black, 1, 0.2))
            scatter!(μ[4:2:end],μ[5:2:end],
            markersize = 1 )
            Plots.gui()
            h5open(include_dir * "iter_$i.h5", "w") do file
                write(file, "/poses", poses)
                write(file, "/μ", μ)
                write(file, "/ξ", ξ)
                write(file, "/Ω", Ω)
                write(file, "/m₀", m₀)
                #write(file, "/ξ_hist", ξ_hist)
            end
            poses = zeros(3,5000)
            #ξ_hist = zeros(3,500)
        end

        motion(speed[i],steering[i],dt,m₀,μ,ξ,Ω,graph)

        estimate(m₀,μ,ξ,Ω)

        m₃ = zeros(Int,0)
        while j ≤ llsr && timelsr[j] < timeit[i]+dt*1000
             z₁ = z[2*j-1:2*j,z[2*j,:] .≠ 0]
             if !isempty(z₁)
                 c = correspondence(z₁,m₀,μ,Ω,graph)
                 Ω,graph = measurement(z₁,c,μ,ξ,Ω,graph)
                 m₃ = sort(union(m₃,c))
             end
             j += 1
        end
        if !isempty(m₃)
                m₄ = sort(setdiff(m₀,m₃))
                q = [m₄;unique(m₃)]
                tmp = q[max(1,length(q)-n+1):end]
                m₁ = setdiff(m₀,tmp)
                m₁ = union(m₁,setdiff(m₃,tmp))
                m₀ = tmp
                if !isempty(m₁)
                    sparsification(m₀,m₁,μ,ξ,Ω,graph)
                end
        end

        t₃ = time()
        println("\33[2J")
        @printf "iter: %d\n" i
        @printf "iter time: %.5f\n" t₃-t₂
        @printf "avg: %.5f\n" ((t₃-t₁)/(i-stindex+1))

        poses[:,1+i%5000] = μ[1:3]
    end
end #function seif()
