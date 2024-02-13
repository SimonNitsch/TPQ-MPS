#include <map>
#include <string>
#include <stdexcept>
#pragma once

namespace TPQ_MPS{

class Hamiltonian{

    private:
    std::map<std::string, double> H;

    public:
    Hamiltonian(){
        H["K"] = 0;
        H["J"] = 0;
        H["hx"] = 0;
        H["hy"] = 0;
        H["hz"] = 0;
        H["Gamma"] = 0;
        H["GammaQ"] = 0;
    }

    double get(std::string x){
        if (H.count(x)){
            return H[x];
        }
        else{
            throw std::invalid_argument("Input Argument is not valid");
        }
    }

    void set(std::string x, double n){
        if (H.count(x)){
            H[x] = n;
        }
        else{
            throw std::invalid_argument("Input Argument is not valid");
        }
    }


};


}


