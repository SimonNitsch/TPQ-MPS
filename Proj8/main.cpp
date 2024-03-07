#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TPQ-MPS-v3-05/main.h"
#include <ios>



int main(){
    std::ios_base::sync_with_stdio(false);
    TPQ_MPS::Hamiltonian H_Details;
    H_Details.set("K",4./3.);
    //H_Details.set("hz",0.15);
    //H_Details.set("J",1);

    int LX = 8;
    int LY = 6;
    int auxiliaries = 7;

    itensor::SpinHalf sites;
    std::string filename = "EPaper";

    int TimeSteps = 100;
    int Evols = 30;
    std::vector<double> intervals = {5,20,35,40};
    int init_states = 256;


    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,sites,auxiliaries,"Honeycomb");
    Model.Time_Evolution(TimeSteps,intervals,Evols,init_states);
    Model.Save(filename);



}


