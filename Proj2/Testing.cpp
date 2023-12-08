#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <chrono>
#include <complex>
#include <ios>
#include <type_traits>


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


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,2> Create_Heisenberg_Model_1D(int N, double J, double beta, T& sites, int auxiliaries){
    sites = T(N+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1+auxiliaries; i != N+auxiliaries; i++){
        ampo += J,"Sz",i,"Sz",i+1;
        ampo += J*0.5,"S+",i,"S-",i+1;
        ampo += J*0.5,"S-",i,"S+",i+1;
    }
    auto H = itensor::toMPO(ampo);
    auto U = itensor::toExpH(ampo,beta/2);

    std::array<itensor::MPO,2> output;
    output[0] = H;
    output[1] = U;
    return output;
}



std::array<int,3> get_neighbour_data(int N, int M, int pos){
    int X = (pos-1) % N;
    int Y = (pos-1) / N;
    std::array<int,3> neighbours{};

    if (X != (N-1)){
        neighbours[0] = pos+1;
    }
    if (Y != 0){
        neighbours[1] = pos-N+1;
    }
    if (pos != 1 && pos != M*N){
        neighbours[2] = pos-1;
    }
    return neighbours;
    
}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,2> Create_Kitaev_Honeycomb_Model_2D(int N, int M, double Kx, double Ky, double Kz, double beta, T& sites, int auxiliaries){
    sites = T(N*M+2*auxiliaries, {"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    std::vector<int> full_points;
    full_points.reserve(((N+1)/2)*M);
    for (int m = 0; m != M; m++){
        for (int n = 1; n <= N; n+=2){
            full_points.push_back(n+N*m);
        }   
    }

    for (int& i : full_points){
        std::array<int,3> neighbours = get_neighbour_data(N,M,i);
        if (neighbours[0] != 0){
            ampo += 2*Kx,"Sx",i+auxiliaries,"Sx",neighbours[0]+auxiliaries;
        }
        if (neighbours[1] != 0){
            ampo += 2*Ky,"Sy",i+auxiliaries,"Sy",neighbours[1]+auxiliaries;
        }
        if (neighbours[2] != 0){
            ampo += 2*Kz,"Sz",i+auxiliaries,"Sz",neighbours[2]+auxiliaries;
        }
    }

    auto H = itensor::toMPO(ampo);
    auto U = itensor::toExpH(ampo,beta/2);

    itensor::PrintData(ampo);


    std::array<itensor::MPO,2> output;
    output[0] = H;
    output[1] = U;
    return output;

}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,2> Create_Heisenberg_Model_2D(int N, int M, double Jx, double Jy, double beta, T& sites, int auxiliaries){
    sites = T(N*M+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

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
        std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));

        for (int j = 0; j != TimeSteps; j++){
            psi = itensor::applyMPO(U,psi,{"Method=","DensityMatrix","MaxDim=",512,"Cutoff=",1e-8});
            std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::innerC(psi,psi);
            E_vec.push_back(std::real(E));
            std::cout << E << "\n";
        }
        Energies.push_back(E_vec);
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << ", Time Needed: " << time.count() << " seconds\n";

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
    std::ios_base::sync_with_stdio(true);
    auto start = std::chrono::system_clock::now();
    
    int N = 6;
    int M = 8;
    double J = 1.;
    double beta = 0.01;
    double K = 1./3.;

    int TimeSteps = 1000;
    int Evols = 25;
    
    itensor::SpinHalf sites;
    auto H_U = Create_Kitaev_Honeycomb_Model_2D(N,M,K,K,K,beta,sites,5);

    std::vector<std::vector<double>> Energies = Calculate_Energies(TimeSteps,Evols,H_U,sites);
    std::vector<double> means = Mean(Energies);

    auto end = std::chrono::system_clock::now();
    auto total_time = std::chrono::duration<double>(end-start);
    std::cout << "Time needed for the whole calculation: " << total_time.count() << " seconds\n";

    for (auto& m : means){
        std::cout << m/((N-1)*M+(N/2+1)*(M-1)) << "\n";
    }



}
