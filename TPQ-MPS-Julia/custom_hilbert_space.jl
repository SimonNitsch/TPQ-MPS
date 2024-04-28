using ITensors



function ITensors.space(::SiteType"Nitsch";size=2)
    return size
end

function ITensors.op(::OpName"Sz", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sz = zeros(ComplexF64,dims,dims)

    for i = 1:dims
        cs = i - dims - 1
        Sz[i,i] = cs
    end
    return ITensor(Sz,prime(s),dag(s))
end

function ITensors.op(::OpName"Sx", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sx = zeros(ComplexF64,dims,dims)

    for i = 1:(dims-1)
        cs = i - dims - 1
        f = sqrt(dims * (dims + 1) - cs * (cs + 1)) / 2.
        Sx[i+1,i] = f
        Sx[i,i+1] = f
    end
    return ITensor(Sx,prime(s),dag(s))
end

function ITensors.op(::OpName"Sy", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sy = zeros(ComplexF64,dims,dims)

    for i = 1:(dims-1)
        cs = i - dims - 1
        f = sqrt(dims * (dims + 1) - cs * (cs + 1)) / 2. * -1im
        Sy[i+1,i] = f
        Sy[i,i+1] = -f
    end
    return ITensor(Sy,prime(s),dag(s))
end

function ITensors.op(::OpName"Expz", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sz = zeros(ComplexF64,dims,dims)

    for i = 1:dims
        cs = i - dims - 1
        Sz[i,i] = cs
    end
    Expz = exp(1im * pi * Sz)
    return ITensor(Expz,prime(s),dag(s))
end

function ITensors.op(::OpName"Expx", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sx = zeros(ComplexF64,dims,dims)

    for i = 1:(dims-1)
        cs = i - dims - 1
        f = sqrt(dims * (dims + 1) - cs * (cs + 1)) / 2.
        Sx[i+1,i] = f
        Sx[i,i+1] = f
    end
    Expx = exp(1im * pi * Sx)
    return ITensor(Expx,prime(s),dag(s))
end

function ITensors.op(::OpName"Expy", ::SiteType"Nitsch", s::Index)
    dims = dim(s)
    Sy = zeros(ComplexF64,dims,dims)

    for i = 1:(dims-1)
        cs = i - dims - 1
        f = sqrt(dims * (dims + 1) - cs * (cs + 1)) / 2. * -1im
        Sy[i+1,i] = f
        Sy[i,i+1] = -f
    end
    Expy = exp(1im * pi * Sy)
    return ITensor(Expy,prime(s),dag(s))
end






