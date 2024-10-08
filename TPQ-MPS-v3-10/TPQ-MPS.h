#include <iostream>
#include <array>
#include <vector>
#include <chrono>
#include "itensor/all.h"
#include <ios>
#include <type_traits>
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>
#include <cstdlib>
#include "TDVP/tdvp.h"
#include <fstream>
#include "Hamiltonian.h"
#include <exception>
#include <filesystem>
#include "customspin_nitsch.h"
#pragma once


std::vector<double> operator-(double x, std::vector<double>& v){
    std::vector<double> vx;
    vx.reserve(v.size());
    
    for (auto& i : v){
        vx.push_back(x-i);
    }
    return vx;
}


using namespace itensor;
namespace TPQ_MPS{



class Kitaev_Model{

    public:
    AutoMPO ampo;
    int Lattice_Type; // 1 for Honeycomb, 2 for Triangular, 3 for Periodic Honeycomb
    bool Calc_Type; // true for DMRG, false for TDVP
    bool Calculate_Susceptibility;

    public:
    MPO H0;
    std::array<MPO,3> M;
    std::array<MPO,3> M2;
    Hamiltonian H_Details;
    CustomSpinNitsch sites;
    CustomSpin sitestan;
    MPO H_flux;
    MPS Htan;

    private:
    double GSE;
    std::vector<std::array<double,2>> E;
    std::vector<std::array<double,2>> Cv;
    std::vector<std::array<double,2>> S;
    std::vector<std::array<double,2>> W;
    std::vector<std::array<double,2>> Chix;
    std::vector<std::array<double,2>> Chiy;
    std::vector<std::array<double,2>> Chiz;


    std::array<int,3> get_neighbour_data_hex(int LX, int LY, int pos);
    std::array<int,3> get_neighbour_data_hex_periodic(int LX, int LY, int pos);
    std::array<int,3> get_neighbour_data_tri(int LX, int LY, int pos);

    std::array<int,3> get_neighbour_data(int LX, int LY, int pos){
        std::array<int,3> n;
        if (Lattice_Type == 1){
            n = get_neighbour_data_hex(LX,LY,pos);
        } else if (Lattice_Type == 2){
            n = get_neighbour_data_tri(LX,LY,pos);
        } else if (Lattice_Type == 3){
            n = get_neighbour_data_hex_periodic(LX,LY,pos);
        }
        return n;
    }


    void add_kitaev_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_magnetic_interaction(int LX, int LY, int aux);
    void add_heisenberg_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gamma_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gammaq_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    MPO honeycomb_flux_operator(int LX, int LY, int aux);
    std::array<std::array<MPO,3>,2> magnetization_operators(int LX, int LY, int aux);

    std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M);
    
    void tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, std::array<std::vector<double>,3>& Chi_vec, MPS& psi, Cplx& t, int TimeSteps, Args& args, Sweeps& Sweeps, double& cb);

    void chi_int(MPS& psi, double n, double t, std::array<std::vector<double>,3>& chi_vec, double step, Sweeps& sweeps, Args& args);

    template<std::size_t n, typename T>
    void save_data(std::string filename, std::vector<std::array<T,n>>& vec);
    
    template<typename T>
    void save_data(std::string filename, T v);

    std::vector<double> derivative(std::vector<double>& f, double dx);
    std::vector<double> integral(std::vector<double>& f, double dx, double c);
    std::vector<double> multiply(std::vector<double>& a, std::vector<double>& b);

    int dims;

    public:
    MPS mpo_to_tanmps(MPO& H);
    MPO mpo_to_tanmpo(MPO& H);
    MPS mps_to_tanmps(MPS& X);
    void tan_tdvp_loop(int steps, double dt, MPS& Hexptan, MPO& H0tan, Sweeps& sweeps, Args& tdvp_args, std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, MPO& Hfluxtan, double& cb);
    double tan_energy(MPS& Hexptan, MPO& H0tan, int states=25, int statesdim=64);

    void Tan_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, int max_sites=256);


    public:
    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, int DoubleSpin, int aux, std::string shape){
        this -> sites = CustomSpinNitsch((LX*LY+2*aux),{"2S=",DoubleSpin,"ConserveQNs=",false});
        this -> ampo = AutoMPO(sites);
        this -> H_Details = H_Details;
        this -> dims = DoubleSpin + 1;
        this -> sitestan = CustomSpin((LX*LY+2*aux),{"2S=",(dims*dims-1),"ConserveQNs=",false});

        std::vector<int> full_points;
        if (shape == "Honeycomb" || shape == "HoneycombPeriodic"){
            if (shape == "Honeycomb"){
                Lattice_Type = 1;
            } else {
                Lattice_Type = 3;
            }
            full_points.reserve(((LY+1)/2)*LX);
            for (int m = 0; m != LX; m++){
                for (int n = 1; n <= LY; n+=2){
                    full_points.push_back(n+LY*m);
                }   
            }
        } else if (shape == "Triangular"){
            Lattice_Type = 2;
            full_points.reserve(LX*LY);
            for (int i = 1; i != LX*LY+1; i++){
                full_points.push_back(i);
            }

        } else {
            std::invalid_argument("Invalid Lattice Shape");
        }

        add_kitaev_interaction(LX,LY,full_points,aux);
        add_magnetic_interaction(LX,LY,aux);
        add_heisenberg_interaction(LX,LY,full_points,aux);
        add_gamma_interaction(LX,LY,full_points,aux);
        add_gammaq_interaction(LX,LY,full_points,aux);


        this -> H_flux = honeycomb_flux_operator(LX,LY,aux);

        //PrintData(H_flux);
        
        this -> H0 = toMPO(this -> ampo);
        //auto Ms = magnetization_operators(LX,LY,aux);
        //M = Ms[0];
        //M2 = Ms[1];

        std::cout << "Spin " << DoubleSpin << "/2 System" << "\n";
        std::cout << "Lattice Type: " << shape << "\n";
        std::cout << "LX = " << LX << ", LY = " << LY << "\n";
        std::cout << "Auxiliary Sites: " << aux << "\n\n";
        std::cout << "Hamiltonian Parameters \n";
        this->H_Details.print();
        std::cout << "\n\n";


        

    }


    void Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, bool Calculate_Suceptibility=false, int max_sites=256, int init_rand_sites=32, std::string TDVP_Type="TwoSite");
    std::array<std::vector<std::array<double,2>>,2> Calculate_Heat_Capacity(int TimeSteps, std::vector<double>& intervals, std::vector<std::vector<double>>& Energies);

    Hamiltonian Get_Constants(){
        return H_Details;
    }

    void Print_Interactions(){
        PrintData(ampo);
    }

    void Save(std::string x){
        if (Calc_Type){
            std::string xGSE = x + "_GSE";
            save_data(xGSE,GSE);
        } 
        else {
            std::filesystem::create_directory(x);
            std::string xE = x + "/" + "E";
            std::string xC = x + "/" + "C";
            std::string xS = x + "/" + "S";
            std::string xW = x + "/" + "W";

            save_data(xE,E);
            save_data(xC,Cv);
            save_data(xS,S);
            save_data(xW,W);

            if (Calculate_Susceptibility){
                std::string xcx = x + "/" + "Chix";
                std::string xcy = x + "/" + "Chiy";
                std::string xcz = x + "/" + "Chiz";

                save_data(xcx,Chix);
                save_data(xcy,Chiy);
                save_data(xcz,Chiz);
            }
        }
        
    }

    
    MPS DMRG(int Sweeps=10){
        Calc_Type = true;
        auto psi0 = randomMPS(sites,16);
        auto [energy, psi] = dmrg(H0,psi0,Sweeps,{"Quiet",true});
        GSE = energy;
        return psi;
    }





};






}










