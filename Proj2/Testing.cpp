#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <chrono>
#include <complex>
#include <ios>
#include <type_traits>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <algorithm>


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

std::array<itensor::MPO,3> Create_Heisenberg_Model_1D(int N, double J, double beta, T& sites, int auxiliaries){
    sites = T(N+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1+auxiliaries; i != N+auxiliaries; i++){
        ampo += J,"Sz",i,"Sz",i+1;
        ampo += J*0.5,"S+",i,"S-",i+1;
        ampo += J*0.5,"S-",i,"S+",i+1;
    }
    auto H = itensor::toMPO(ampo);
    auto U1 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
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

std::array<itensor::MPO,3> Create_Kitaev_Honeycomb_Model_2D(int N, int M, double Kx, double Ky, double Kz, double beta, T& sites, int auxiliaries){
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
            ampo += Kx,"Sx",i+auxiliaries,"Sx",neighbours[0]+auxiliaries;
        }
        if (neighbours[1] != 0){
            ampo += Ky,"Sy",i+auxiliaries,"Sy",neighbours[1]+auxiliaries;
        }
        if (neighbours[2] != 0){
            ampo += Kz,"Sz",i+auxiliaries,"Sz",neighbours[2]+auxiliaries;
        }
    }

    itensor::PrintData(ampo);
    auto H = itensor::toMPO(ampo);
    auto U1 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
    return output;

}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,3> Create_Heisenberg_Model_2D(int N, int M, double J, double beta, T& sites, int auxiliaries){
    sites = T(N*M+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != N; i++){
        for (int j = 0; j != M; j++){
            int n = auxiliaries+j*N+i;
            ampo += J,"Sz",n,"Sz",n+1;
            ampo += J*0.5,"S+",n,"S-",n+1;
            ampo += J*0.5,"S-",n,"S+",n+1;
        }
    }

    for (int i = 1; i <= N; i++){
        for (int j = 0; j != (M-1); j++){
            int n = auxiliaries+j*N+i;
            ampo += J,"Sz",n,"Sz",n+N;
            ampo += J*0.5,"S+",n,"S-",n+N;
            ampo += J*0.5,"S-",n,"S+",n+N;
        }
    }
    
    itensor::PrintData(ampo);
    auto H = itensor::toMPO(ampo);
    auto U1 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
    return output;


}





std::vector<std::vector<double>> Calculate_Energies(int TimeSteps, int Evols, std::array<itensor::MPO,3>& H_U, itensor::SiteSet& sites){
    auto H = H_U[0];
    auto U1 = H_U[1];
    auto U2 = H_U[2];
    
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);
    std::vector<double> E_vec;
    E_vec.reserve(TimeSteps+1);
    auto IState = itensor::InitState(sites);


    for (int i = 0; i != Evols; i++){
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(IState);
        //std::cout << "Hello\n";
        std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::inner(psi,psi);
        //std::cout << E << "\n";
        //std::cout << Ereal << "\n";
        E_vec.push_back(std::real(E));
        //std::cout << "b" << std::flush;

        for (int j = 0; j != TimeSteps; j++){
            //std::cout << "a";
            psi = itensor::applyMPO(U1,psi,{"Method=","DensityMatrix","MaxDim=",256,"Cutoff=",1e-8});
            psi = itensor::applyMPO(U2,psi,{"Method=","DensityMatrix","MaxDim=",256,"Cutoff=",1e-8});
            std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::innerC(psi,psi);
            E_vec.push_back(std::real(E));
            //std::cout << "a";
            //std::cout << j << std::flush;
            //std::cout << E << "\n";
        }
        Energies.push_back(E_vec);
        E_vec.clear();
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << ", Time Needed: " << time.count() << " seconds\n" << std::flush;

    }

    return Energies;

}


std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M){
    int M0 = M[0].size();
    int M1 = M.size();
    std::vector<std::array<double,2>> vec;
    vec.reserve(M0);
    
    for (int i = 0; i != M0; i++){
        double v = 0;
        double v2 = 0;
        for (int j = 0; j != M1; j++){
            v += M[j][i];
            v2 += M[j][i] * M[j][i];
        }
        double mean = v/static_cast<double>(M1);
        double std = std::sqrt(v2/M1 - mean*mean);
        std::array<double,2> v = {mean,std};
        vec.emplace_back(v);
    }
    
    return vec;

}


template<std::size_t n, typename T>
void Save_Data(std::string& filename, std::vector<std::array<T,n>>& vec, int data_points=100){
    std::ofstream file(filename,std::ios::binary);
    int length = std::min(vec.size(),data_points);

    for (int i = 0; i != length; i++){
        int next_index = (i*vec.size())/length;
        std::array<T,n> point = vec[next_index];
        file.write(reinterpret_cast<const char*>(point.data()),point.size()*sizeof(T));
    }
    file.close();

    std::cout << "Data saved as: " << filename << "\n" << std::flush;

}




int main(){
    std::ios_base::sync_with_stdio(false);
    auto start = std::chrono::system_clock::now();
    
    int N = 6;
    int M = 8;
    double J = 1.;
    double beta = 0.05;
    double K = 1.;

    int TimeSteps = 300;
    int Evols = 1;
    
    itensor::SpinHalf sites;
    auto H_U = Create_Kitaev_Honeycomb_Model_2D(N,M,K,K,K,beta,sites,7);

    std::vector<std::vector<double>> Energies = Calculate_Energies(TimeSteps,Evols,H_U,sites);
    std::vector<std::array<double,2>> means = Mean(Energies);

    auto end = std::chrono::system_clock::now();
    auto total_time = std::chrono::duration<double>(end-start);
    std::cout << "Time needed for the whole calculation: " << total_time.count() << " seconds\n";

    for (auto& i : means){
        std::cout << i[0] << " " << i[1] << "\n";
    }



}
