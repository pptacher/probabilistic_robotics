using Distributed
const N_CORES = 4
rmprocs(workers())
addprocs(N_CORES)
@everywhere push!(LOAD_PATH,".")
@everywhere using SparseArrays, LinearAlgebra, SharedArrays
include("seif.jl")

seif()
