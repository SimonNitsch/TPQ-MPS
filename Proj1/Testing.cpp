#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <chrono>
#include <complex>
#include <ios>


std::array<itensor::MPO,2> Create_Heisenberg_Model_1D(int N, double J, double beta, itensor::SpinHalf& sites){
    sites = itensor::SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != N; i++){
        ampo += J,"Sz",i,"Sz",i+1;
        ampo += J*0.5,"S+",i,"S-",i+1;
        ampo += J*0.5,"S-",i,"S+",i+1;
    }
    auto H = itensor::toMPO(ampo);
    auto U = itensor::toExpH(ampo,beta);

    std::array<itensor::MPO,2> output;
    output[0] = H;
    output[1] = U;
    return output;
}

std::vector<std::vector<double>> Calculate_Energies(int TimeSteps, int Evols, std::array<itensor::MPO,2>& H_U, itensor::SiteSet& sites){
    auto H = H_U[0];
    auto U = H_U[1];
    
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);
    std::vector<double> E_vec;
    E_vec.reserve(TimeSteps+1);
    auto IState = itensor::InitState(sites);

    for (int i = 0; i != Evols; i++){
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(IState);
        std::complex<double> E = itensor::inner(psi,H,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));

        for (int j = 0; j != TimeSteps; j++){
            psi = itensor::applyMPO(U,psi,{"Method=","DensityMatrix","MaxDim=",20,"Cutoff=",1e-5});
            std::complex<double> E = itensor::inner(psi,H,psi) / itensor::inner(psi,psi);
            E_vec.push_back(std::real(E));
        }
        Energies.push_back(E_vec);
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << ", Time Needed: " << time.count() << "seconds\n";

    }

    return Energies;

}

std::vector<double> Mean(std::vector<std::vector<double>>& M){
    std::vector<double> vec;
    int M0 = M[0].size();
    int M1 = M.size();
    vec.reserve(M0);
    
    for (int i = 0; i != M0; i++){
        double v = 0;
        for (int j = 0; j != M1; j++){
            v += M[j][i];
        }
        vec.push_back(v/M1);

    }
    return vec;

}




int main(){
    std::ios_base::sync_with_stdio(false);

    int N = 10;
    double J = 1.;
    double beta = 0.05;

    int TimeSteps = 150;
    int Evols = 25;

    itensor::SpinHalf sites;
    auto H_U = Create_Heisenberg_Model_1D(10,1.,0.05,sites); 

    std::vector<std::vector<double>> Energies = Calculate_Energies(TimeSteps,Evols,H_U,sites);
    std::vector<double> means = Mean(Energies);

    for (auto& m : means){
        std::cout << m << "\n";
    }



}
