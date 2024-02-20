#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TPQ-MPS-v3/TPQ-MPS.h"
#include "TPQ-MPS-v3/TPQ-MPS.cpp"
#include <ios>



int main(){
    //std::ios_base::sync_with_stdio(false);
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
    std::vector<double> intervals = {1};

    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,sites,auxiliaries,"Honeycomb");
    Model.Print_Interactions();
    Model.Time_Evolution(TimeSteps,intervals,Evols);
    Model.Save(filename);


}




