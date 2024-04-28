using ITensors
include("custom_hilbert_space.jl")





function create_mpo(LX::Int, LY::Int, aux::Int, size::Int, H_Details::Dict, Lattice::String, posarr::Array)

    os = OpSum()
    os = add_kitaev_interaction!(LX,LY,aux,H_Details,Lattice,os,posarr)
    os = add_magnetic_interaction!(LX,LY,aux,H_Details,os)
    os = add_heisenberg_interaction!(LX,LY,aux,H_Details,Lattice,os,posarr)
    os = add_gamma_interaction!(LX,LY,aux,H_Details,Lattice,os,posarr)
    os = add_gammaq_interaction!(LX,LY,aux,H_Details,Lattice,os,posarr)
    l::Int = LX * LY + 2 * aux
    @show os

    sites = siteinds("Nitsch",l; size)
    H = MPO(ComplexF64,os,sites)
    return (H,sites)
end


function create_flux_mpo(LX::Int, LY::Int, aux::Int, sites::Vector{<:Index})
    PY::Int = LY/2 - 1
    PX::Int = LX - 1
    Wfac = 64. / convert(Float64,PY) / convert(Float64,PX)
    fos = OpSum()

    for i = 1:PX
        for j = 1:PY
            f = LY*(i-1) + 2*j + aux
            fos += Wfac,"Expx",f,"Expy",f+1,"Expz",f+2,"Expx",f+LY+1,"Expy",f+LY,"Expz",f+LY-1
            println("a")
        end
    println("b")
    end

    H_flux = MPO(ComplexF64,fos,sites)
    println("c")
    return H_flux
end
