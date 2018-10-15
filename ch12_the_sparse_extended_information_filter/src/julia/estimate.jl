function estimate(m₀,μ,ξ,Ω)
    n = (length(μ)-3)÷2

    if n < 100
        μ .= Ω \ ξ
    else
        @inbounds for i ∈ m₀
            ind = [2*i+2,2*i+3]
            μ[ind] = Ω[ind,ind] \ (ξ[ind]-Ω[ind,:]*μ+Ω[ind,ind]*μ[ind])
        end
        μ[1:3] = Ω[1:3,1:3] \ (ξ[1:3]-Ω[1:3,:]*μ+Ω[1:3,1:3]*μ[1:3])
    end
end
