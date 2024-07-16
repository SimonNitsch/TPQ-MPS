#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TPQ-MPS-v3-13/main.h"
#include <ios>
#include <stdlib.h>



int main(){
    char env[] = "MKL_NUM_THREADS=16";
    putenv(env);
    std::ios_base::sync_with_stdio(false);
    TPQ_MPS::Hamiltonian H_Details;
    H_Details.set("K",1.);
    //H_Details.set("K",0.809017);
    //H_Details.set("J",0.587785);
    //H_Details.set("hz",0.005);
    //H_Details.set("J",1);

    int LX = 2;
    int LY = 8;
    int auxiliaries = 0;

    int S2 = 1;
    std::string filename = "2by8";

    int max_sites = 2048;
    int min_sites = 32;

    int TimeSteps = 70;
    std::vector<double> intervals = {5,8,12,15,20};

    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,S2,auxiliaries,"HoneycombPeriodic");
    Model.Tan_Evolution(TimeSteps,intervals,max_sites);
    //auto psi = Model.DMRG();
    Model.Save(filename);



}


