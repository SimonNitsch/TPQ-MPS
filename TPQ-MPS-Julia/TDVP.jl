using ITensors
using ITensorTDVP
using Statistics

struct MPO_Object
    H::MPO
    H_flux::MPO
    sites::Vector{<:Index}
    length::Int
end


function tdvp_loop(X::MPO_Object, psi::MPS, interval::Float64, step::Int, EMat::Matrix, CMat::Matrix, SMat::Matrix, FMat::Matrix, curind::Int, curevol::Int, curbeta::Float64, maxdim::Int)
    tstep = -interval / convert(Float64,step) / 2.
    
    for i = 1:step
        curbeta -= tstep * 2
        psi = tdvp(X.H,-tstep,psi;nsweeps=1,normalize=true,maxdim=maxdim)

        EMatC::ComplexF64 = inner(psi',X.H,psi)
        EMat[curind+i,curevol] = real(EMatC)
        CMat[curind+i,curevol] = (EMat[curind+i,curevol] - EMat[curind+i-1,curevol]) / (-2 * tstep) * curbeta * curbeta
        SMat[curind+i,curevol] = SMat[curind+i-1,curevol] + (CMat[curind+i,curevol] / (curbeta + 2 * tstep) + CMat[curind+i-1,curevol] / curbeta) * -tstep
        FMatC::ComplexF64 = inner(psi',X.H_flux,psi)
        FMat[curind+i,curevol] = real(FMatC)
        println("Time Step Done")

    end
    return (psi,EMat,CMat,SMat,FMat,curbeta)
end


function matmean(x::Matrix)
    xs1 = size(x,1)
    meanx = zeros(xs1,2)
    
    for i = 1:xs1
        xvec = x[i,:]
        meanx[i,1] = mean(xvec)
        meanx[i,2] = std(xvec)
    end
    return meanx
end






