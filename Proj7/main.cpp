#include<iostream>
#include<array>
#include<vector>
#include "itensor/all.h"
#include "TPQ-MPS-v2-3/TPQ-MPS.h"




int main(){
    int N = 8;
    int M = 6;
    itensor::SpinHalf sites;

    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {0,0,0.15};
    double J = 1;
    int auxiliaries = 7;

    std::vector<int> intervals = {5,10,20};
    int TimeSteps = 120;
    int Evols = 50;

    auto H = TPQ_MPS::Create_Kitaev_Honeycomb_Model_2D(N,M,K,J,h,sites,auxiliaries);

}






