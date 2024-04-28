using ITensors


function get_neighbour_data_hex(pos::Int, LX::Int, LY::Int)
    X = (pos - 1) / LX
    Y = (pos - 1) % LY
    neighbours = zeros(Int,3)

    if Y != (LY - 1)
        neighbours[1] = pos + 1
    end
    if X != 0
        neighbours[2] = pos - LY + 1
    end
    if pos != 1 && pos != LX*LY
        neighbours[3] = pos - 1
    end
    return neighbours
end

function get_neighbour_data_hex_per(pos::Int, LX::Int, LY::Int)
    X = (pos - 1) / LX
    Y = (pos - 1) % LY
    neighbours = zeros(Int,3)

    if Y != (LY - 1)
        neighbours[1] = pos + 1
    end
    if X != 0
        neighbours[2] = pos - LY + 1
    else
        neighbours[2] = pos + (LX-1) * LY + 1
    end
    if Y != 0
        neighbours[3] = pos - 1
    else
        neighbours[3] = pos + LY - 1
    end
    return neighbours
end

function get_neighbour_data_tri(pos::Int, LX::Int, LY::Int)
    X = (pos-1) / LY;
    Y = (pos-1) % LY;
    neighbours = zeros(Int,3);

    if X != LX - 1
        neighbours[1] = pos + LY
        if Y == 0
            neighbours[3] = pos + 2 * LY - 1
        else
            neighbours[3] = pos + LY - 1
        end
    else
        neighbours[1] = pos - (LX-1) * LY
        if Y == 0
            neighbours[3] = pos - (LX-2) * LY - 1
        else
            neighbours[3] = pos - (LX-1) - 1
        end
    end

    if Y == LY - 1
        neighbours[2] = pos - LY + 1
    else
        neighbours[2] = pos + 1
    end
    return neighbours    
end


function  get_neighbour_data(pos::Int, LX::Int, LY::Int, Lattice::String)
    if Lattice == "Honeycomb"
        neighbours = get_neighbour_data_hex(pos,LX,LY)
    elseif Lattice == "HoneycombPeriodic"
        neighbours = get_neighbour_data_hex_per(pos,LX,LY)
    elseif Lattice == "Triangular"
        neighbours = get_neighbour_data_tri(pos,LX,LY)
    end
    return neighbours    
end


function add_kitaev_interaction!(LX::Int, LY::Int, aux::Int, H::Dict, Lattice::String, opsum::OpSum, posarr::Array)
    l = length(posarr)
    @show posarr
    Kx = getindex(H,"Kx")
    Ky = getindex(H,"Ky")
    Kz = getindex(H,"Kz")

    for i = 1:l
        n = get_neighbour_data(posarr[i],LX,LY,Lattice)
        if n[1] != 0
            opsum += Kx,"Sx",posarr[i]+aux,"Sx",n[1]+aux
        end
        if n[2] != 0
            opsum += Ky,"Sy",posarr[i]+aux,"Sy",n[2]+aux
        end
        if n[3] != 0
            opsum += Kz,"Sz",posarr[i]+aux,"Sz",n[3]+aux
        end
    end
    return opsum
end


function add_heisenberg_interaction!(LX::Int, LY::Int, aux::Int, H::Dict, Lattice::String, opsum::OpSum, posarr::Array)
    l = length(posarr)
    J = getindex(H,"J")

    for i = 1:l 
        n = get_neighbour_data(posarr[i],LX,LY,Lattice)
        if n[1] != 0
            opsum += J,"Sx",posarr[i]+aux,"Sx",n[1]+aux
            opsum += J,"Sy",posarr[i]+aux,"Sy",n[1]+aux
            opsum += J,"Sz",posarr[i]+aux,"Sz",n[1]+aux
        end
        if n[2] != 0
            opsum += J,"Sx",posarr[i]+aux,"Sx",n[2]+aux
            opsum += J,"Sy",posarr[i]+aux,"Sy",n[2]+aux
            opsum += J,"Sz",posarr[i]+aux,"Sz",n[2]+aux
        end
        if n[3] != 0
            opsum += J,"Sx",posarr[i]+aux,"Sx",n[3]+aux
            opsum += J,"Sy",posarr[i]+aux,"Sy",n[3]+aux
            opsum += J,"Sz",posarr[i]+aux,"Sz",n[3]+aux
        end        
    end
    return opsum
end


function add_magnetic_interaction!(LX::Int, LY::Int, aux::Int, H::Dict, opsum::OpSum)
    hx = getindex(H,"hx")
    hy = getindex(H,"hy")
    hz = getindex(H,"hz")

    for i = 1:(LX*LY)
        opsum += hx,"Sx",i+aux
        opsum += hy,"Sy",i+aux
        opsum += hz,"Sz",i+aux
    end
    return opsum
end


function add_gamma_interaction!(LX::Int, LY::Int, aux::Int, H::Dict, Lattice::String, opsum::OpSum, posarr::Array)
    l = length(posarr)
    G = getindex(H,"Gamma")

    for i = 1:l
        n = get_neighbour_data(posarr[i],LX,LY,Lattice)
        if n[1] != 0
            opsum += G,"Sy",posarr[i]+aux,"Sz",n[1]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sy",n[1]+aux
        end
        if n[2] != 0
            opsum += G,"Sx",posarr[i]+aux,"Sz",n[2]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sx",n[2]+aux
        end
        if n[3] != 0
            opsum += G,"Sx",posarr[i]+aux,"Sy",n[3]+aux
            opsum += G,"Sy",posarr[i]+aux,"Sx",n[3]+aux
        end
    end
    return opsum
end



function add_gammaq_interaction!(LX::Int, LY::Int, aux::Int, H::Dict, Lattice::String, opsum::OpSum, posarr::Array)
    l = length(posarr)
    G = getindex(H,"Gamma")

    for i = 1:l
        n = get_neighbour_data(posarr[i],LX,LY,Lattice)
        if n[1] != 0
            opsum += G,"Sx",posarr[i]+aux,"Sz",n[1]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sx",n[1]+aux
            opsum += G,"Sx",posarr[i]+aux,"Sy",n[1]+aux
            opsum += G,"Sy",posarr[i]+aux,"Sx",n[1]+aux
        end
        if n[2] != 0
            opsum += G,"Sy",posarr[i]+aux,"Sz",n[2]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sy",n[2]+aux
            opsum += G,"Sy",posarr[i]+aux,"Sx",n[2]+aux
            opsum += G,"Sx",posarr[i]+aux,"Sy",n[2]+aux
        end
        if n[3] != 0
            opsum += G,"Sx",posarr[i]+aux,"Sz",n[3]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sx",n[3]+aux
            opsum += G,"Sy",posarr[i]+aux,"Sz",n[3]+aux
            opsum += G,"Sz",posarr[i]+aux,"Sy",n[3]+aux
        end
    end
    return opsum
end




