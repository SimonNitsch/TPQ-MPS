#include<iostream>
#include<array>
#include<vector>
#include<map>
#include<string>
#include "itensor/all.h"
#include "TPQ-MPS-v3/TPQ-MPS.h"
#include "TPQ-MPS-v3/TPQ-MPS.cpp"




int main(){
    int N = 3;
    int M = 4;
    itensor::SpinHalf sites;

    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {0,0,0.15};
    double J = 1;

    TPQ_MPS::Hamiltonian H_Details;
    H_Details.set("K",1);
    H_Details.set("hz",0.15);
    H_Details.set("J",1);
    int auxiliaries = 7;

    std::vector<double> intervals = {5,10,20};
    int TimeSteps = 120;
    int Evols = 50;

    TPQ_MPS::Kitaev_Model Model(N,M,H_Details,sites,auxiliaries,"Honeycomb");
    Model.Print_Interactions();
    Model.Time_Evolution(TimeSteps,intervals,Evols);
    

}






