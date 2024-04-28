using ITensors
include("opsum_funcs.jl")
include("custom_hilbert_space.jl")
include("TDVP.jl")
include("create_mpo.jl")
include("hamiltonian.jl")





function Create_Kitaev_MPO(LX::Int, LY::Int, auxiliaries::Int, DoubleSpin::Int, H_Details::Dict, Lattice::String)
    length::Int = LX * LY + 2 * auxiliaries
    dims::Int = DoubleSpin + 1

    sitelist = zeros(Int,0)
    if Lattice == "Honeycomb" || Lattice == "HoneycombPeriodic"
        for i = 1:LX
            for j = 1:(LY/2)
                append!(sitelist, (i-1) * LY + j * 2 - 1)
            end
        end
    elseif Lattice == "Triangular"
        for i = 1:(LX*LY)
            append!(sitelist,i)
        end
    end

    H, sites = create_mpo(LX,LY,auxiliaries,dims,H_Details,Lattice,sitelist)
    H_flux = create_flux_mpo(LX,LY,auxiliaries,sites)
    H_Object = MPO_Object(H,H_flux,sites,length)
    return H_Object
end




function TimeEvolution(X::MPO_Object, Intervals::Vector{Float64}, Steps::Vector{Int}, Evolutions::Int; maxdims::Int=256, mindims::Int=32)
    total_steps = sum(Steps)
    EnergyMat = zeros(Float64,total_steps+1,Evolutions)
    CapacityMat = zeros(Float64,total_steps+1,Evolutions)
    EntropyMat = zeros(Float64,total_steps+1,Evolutions)
    FluxMat = zeros(Float64,total_steps+1,Evolutions)
    println("Starting TDVP Algorithm")

    for i = 1:Evolutions
        t1 = time()
        psi = randomMPS(X.sites;linkdims=mindims)
        beta::Float64 = 0.
        curind = 1

        for (in,st) in zip(Intervals,Steps)
            psi, EnergyMat, CapacityMat, EntropyMat, FluxMat, beta = tdvp_loop(X,psi,in,st,EnergyMat,CapacityMat,EntropyMat,FluxMat,curind,i,beta,maxdims)
            curind += st
        end
        elt = time() - t1
        println("Finished Evolution Number: ", i, ", Elapsed Time: ", elt, " s")
    end
    
    Energy = matmean(EnergyMat)
    Capacity = matmean(CapacityMat)
    Entropy = matmean(EntropyMat)
    Flux = matmean(FluxMat)
    return (Energy,Capacity,Entropy,Flux)
end






#=
LX = 3
LY = 4
auxiliaries = 5
DoubleSpin = 2

Intervals = [1.,2.]
Steps = [5,5]
Evolutions = 5

vektor = range(1,5,5)
vektor = 2 .^ vektor
@show vektor

H_Details = Hamiltonian_Variables()
SetH!(H_Details,1.,"K")
Lattice = "HoneycombPeriodic"

M = [0. 1. 0. 0. 
     1. 0. 2. 0.
     0. 2. 0. 1.
     0. 0. 1. 0.]

M2 = exp(pi * 1im * M)
#M2[M2 .< 1e-5] = 0

@show M
@show M2



H = Create_Kitaev_MPO(LX,LY,auxiliaries,DoubleSpin,H_Details,Lattice)
@show H.H
E, C, S, F = TimeEvolution(H,Intervals,Steps,Evolutions)

@show E
@show C
@show S
@show F

=#












