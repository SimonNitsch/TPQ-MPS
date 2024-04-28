


function Hamiltonian_Variables()
    x = Dict("Kx" => 0, "Ky" => 0, "Kz" => 0,
    "J" => 0, "hx" => 0, "hy" => 0, "hz" => 0,
    "Gamma" => 0, "Gamma" => 0)
    return x
end


function SetH!(x::Dict, value::Float64, variable::String)
    if variable=="K"
        setindex!(x,value,"Kx")
        setindex!(x,value,"Ky")
        setindex!(x,value,"Kz")
    else
        setindex!(x,value,variable)
    end

    return x
end

function PrintH!(x::Dict)
    println("Hamiltonian Variables: ")
    println(string("Kx: ",getindex(x,"Kx")))
    println(string("Ky: ",getindex(x,"Ky")))
    println(string("Kz: ",getindex(x,"Kz")))
    println(string("J: ",getindex(x,"J")))
    println(string("hx: ",getindex(x,"hx")))
    println(string("hy: ",getindex(x,"hy")))
    println(string("hz: ",getindex(x,"hz")))
    println(string("Gamma: ", getindex(x,"Gamma")))
    println(string("GammaQ: ", getindex(x,"GammaQ")))
end






