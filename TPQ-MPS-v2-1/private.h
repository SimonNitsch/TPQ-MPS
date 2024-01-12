#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <cmath>
#include "TDVP/tdvp.h"
#include <type_traits>
#pragma once


namespace priv{

void tdvp_loop(std::vector<double>& E_vec, itensor::MPO& H, itensor::MPS& psi, itensor::Cplx t, int Sweeps, int TimeSteps, int data){
    
    int count = 0;
    for (int j = 0; j != TimeSteps; j++){
        double E = itensor::tdvp(psi,H,t,Sweeps,{"Truncate",true,
                                                "DoNormalize",true,
                                                "Quiet",true,
                                                "NumCenter",1,
                                                "ErrGoal",1E-7});
            
        if (j*data >= count*TimeSteps){
            E_vec.push_back(E);
            count++;
        }
    }
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
        std::array<double,2> proto_vec = {mean,std};
        vec.emplace_back(proto_vec);
    }
    
    return vec;

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

itensor::AutoMPO kitaev_honeycomb_2D(int N, int M, std::array<double,3>& K, T& sites, int aux){
    sites = T(N*M+2*aux, {"ConserveQNs=",false});
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
            ampo += K[0],"Sx",i+aux,"Sx",neighbours[0]+aux;
        }
        if (neighbours[1] != 0){
            ampo += K[1],"Sy",i+aux,"Sy",neighbours[1]+aux;
        }
        if (neighbours[2] != 0){
            ampo += K[2],"Sz",i+aux,"Sz",neighbours[2]+aux;
        }
    }

    return ampo;
}


std::array<int,3> get_neighbour_data_tri(int N, int M, int pos){
    int X = (pos-1) / N;
    int Y = (pos-1) % N;
    std::array<int,3> neighbours = {0,0,0};

    if (X != (M-1)){
        neighbours[0] = pos + N;
    }
    if (Y == (N-1)){
        neighbours[1] = pos + 1 - N;
    }
    else{
        neighbours[1] = pos + 1;
    }
    if (X != (M-1)){
        if (Y == 0){
            neighbours[2] = pos + 2*N - 1;
        }
        else{
            neighbours[2] = pos + N - 1;
        }
    }
    return neighbours;

}



template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

itensor::AutoMPO kitaev_tri_2d(int N, int M, std::array<double,3>& K, T& sites, int aux){
    sites = T(N*M+2*aux,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != (N*M+1); i++){
        std::array<int,3> neighbours = get_neighbour_data_tri(N,M,i);
        if (neighbours[0] != 0){
            ampo += K[0],"Sx",i+aux,"Sx",neighbours[0]+aux;
        }
        if (neighbours[1] != 0){
            ampo += K[1],"Sy",i+aux,"Sy",neighbours[1]+aux;
        }
        if (neighbours[2] != 0){
            ampo += K[2],"Sz",i+aux,"Sz",neighbours[2]+aux;
        }
    }

    return ampo; 
    
}




void add_magnetic_field(itensor::AutoMPO& ampo, std::array<double,3>& h, int length, int aux){
    for (int i = (aux+1); i != (length-aux+1); i++){
        ampo += -1*h[0],"Sx",i;
        ampo += -1*h[1],"Sy",i;
        ampo += -1*h[2],"Sz",i;
    }
}







}
