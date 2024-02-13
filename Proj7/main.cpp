#include<iostream>
#include<array>
#include<vector>
#include<map>
#include<string>
#include "itensor/all.h"
#include "TPQ-MPS-v2-3/TPQ-MPS.h"
#include "TPQ-MPS-v3/TPQ-MPS.h"




int main(){
    int N = 4;
    int M = 3;
    itensor::SpinHalf sites;

    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {0,0,0.15};
    double J = 1;
    std::map<std::string,double> H_Details;
    H_Details["K"] = 1;
    H_Details["hz"] = 0.15;
    H_Details["J"] = 1;
    int auxiliaries = 7;

    std::vector<int> intervals = {5,10,20};
    int TimeSteps = 120;
    int Evols = 50;

    auto H = TPQ_MPS::Create_Kitaev_Honeycomb_Model_2D(N,M,H_Details,sites,auxiliaries);

}






