using ITensors
using LinearAlgebra
include("main.jl")
include("custom_hilbert_space.jl")
cpus = Sys.CPU_THREADS
@show cpus
used_cpus = ITensors.blas_get_num_threads()
@show used_cpus
BLAS.set_num_threads(2)
used_cpus = ITensors.blas_get_num_threads()
@show used_cpus



LX = 3
LY = 4
auxiliaries = 5
DoubleSpin = 2

Intervals = [1.,2.]
Steps = [5,5]
Evolutions = 5

H_Details = Hamiltonian_Variables()
SetH!(H_Details,1.,"K")
Lattice = "HoneycombPeriodic"

ampo = OpSum()
ampo += "Sx",1

sites = siteinds("Nitsch",1)
mpo = MPO(ComplexF64,ampo,sites)
@show mpo



#H = Create_Kitaev_MPO(LX,LY,auxiliaries,DoubleSpin,H_Details,Lattice)
#@show H.H
#E, C, S, F = TimeEvolution(H,Intervals,Steps,Evolutions)


