#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TPQ-MPS-v3-01/main.h"
#include <ios>



int main(){
    std::ios_base::sync_with_stdio(false);
    TPQ_MPS::Hamiltonian H_Details;
    H_Details.set("K",1);
    H_Details.set("hx",0.15);
    H_Details.set("J",1);

    int LX = 4;
    int LY = 5;
    int auxiliaries = 7;

    itensor::SpinHalf sites;
    std::string filename = "E";

    int TimeSteps = 20;
    int Evols = 20;
    std::vector<double> intervals = {1,2};

    std::vector<double> u = {0,1,4,16,25,75,0,36,40,20,5,1};
    std::vector<std::vector<double>> U;
    U.reserve(10);
    for (int i = 0; i != 10; i++){
        U.push_back(u);
    }
    int TT = u.size() / intervals.size() - 1;
    std::cout << TT << "\n\n";

    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,sites,auxiliaries,"Honeycomb");
    Model.Print_Interactions();
    //Model.Time_Evolution(TimeSteps,intervals,Evols);
    auto v = Model.Calculate_Heat_Capacity(TT,intervals,U);
    //Model.Save(filename);
    
    for (auto& i : v){
        for (auto& j : i){
            for (auto& k : j){
                std::cout << k << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }


}





