#include<iostream>
#include<array>
#include<vector>
#include<map>
#include<string>
#include "itensor/all.h"
#include "TPQ-MPS-v3/TPQ-MPS.h"
#include "TPQ-MPS-v3/TPQ-MPS.cpp"
#include <ios>




int main(){
    std::ios_base::sync_with_stdio(false);
    int N = 3;
    int M = 4;
    itensor::SpinHalf sites;

    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {0,0,0.15};
    double J = 1;

    TPQ_MPS::Hamiltonian H_Details_H;
    H_Details_H.set("K",1);
    H_Details_H.set("hz",0.15);
    H_Details_H.set("J",1);

    std::map<std::string, double> H_Details_map;
    H_Details_map["K"]  = 1;
    H_Details_map["hz"] = 0.15;
    H_Details_map["J"] = 1;
    int auxiliaries = 7;

    std::vector<double> intervals = {1};
    int TimeSteps = 10;
    int Evols = 1;

    auto Model = TPQ_MPS::Kitaev_Model(M,N,H_Details_H,sites,auxiliaries,"Honeycomb");
    Model.Print_Interactions();
    Model.Time_Evolution(TimeSteps,intervals,Evols); 

    /*
    auto H = TPQ_MPS::Create_Kitaev_Honeycomb_Model_2D(N,M,H_Details_map,sites,auxiliaries);
    auto E = TPQ_MPS::Calculate_Energies_TDVP(TimeSteps,intervals,Evols,H,sites);
    */

}






