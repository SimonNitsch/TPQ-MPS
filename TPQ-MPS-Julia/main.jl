using ITensors
using Plots
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




function TDVP_Time_Evolution(X::MPO_Object, Intervals::Vector{Float64}, Steps::Vector{Int}, Evolutions::Int; maxdims::Int=256, mindims::Int=32)
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






function Plot(E::Matrix{Float64}, C::Matrix{Float64}, S::Matrix{Float64}, F::Matrix{Float64}, Steps::Vector{Int}, Intervals::Vector{Float64})
    x = zeros(Float64,0)
    intv = 0.
    for (in, st) in zip(Intervals,Steps)
        xp = range(intv,length=st+1,stop=intv+in)
        xp = xp[2:end]
        x = [x;xp]
        intv += in
    end
    x = 1. ./ x
    pyplot()

    Em = E[2:end,1]
    Es = E[2:end,2]
    plt1 = plot(x,Em,label="Energy",linecolor=:blue)
    plot!(plt1,x,[Em+Es,Em-Es],linestyle=[:dash],label="",linecolor=:blue)
    plot!(plt1,size=(1920,600))
    plot!(plt1,xaxis=:log)
    xlabel!(plt1,"T")
    ylabel!(plt1,"E")
    png(plt1,"E.png")


    Cm = C[2:end,1]
    Cs = C[2:end,2]
    Sm = S[2:end,1]
    Ss = S[2:end,2]
    plt2 = plot(x,Cm,label="Heat Capacity",linecolor=:red)
    plot!(plt2,x,[Cm+Cs,Cm-Cs],linestyle=[:dash],label="",linecolor=:red)
    plot!(plt2,x,Sm,label="Entropy",linecolor=:green)
    plot!(plt2,x,[Sm+Ss,Sm-Ss],linestyle=[:dash],label="",linecolor=:green)
    plot!(plt2,size=(1920,600))
    plot!(plt2,xaxis=:log)
    xlabel!(plt2,"T")
    ylabel!(plt2,"E")
    png(plt2,"CS.png")


    
    Fm = F[2:end,1]
    Fs = F[2:end,2]
    plt3 = plot(x,Fm,label="Energy",linecolor=:yellow)
    plot!(plt3,x,[Fm+Fs,Fm-Fs],linestyle=[:dash],label="",linecolor=:yellow)
    plot!(plt3,size=(1920,600))
    plot!(plt3,xaxis=:log)
    xlabel!(plt3,"T")
    ylabel!(plt3,"E")
    png(plt3,"F.png")


    return nothing
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












