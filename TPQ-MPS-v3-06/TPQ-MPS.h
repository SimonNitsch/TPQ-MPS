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
#pragma once


namespace TPQ_MPS{



class Kitaev_Model{

    private:
    itensor::AutoMPO ampo;
    int Lattice_Type; // 1 for Honeycomb, 2 for Triangular
    bool Calc_Type; // true for DMRG, false for TDVP

    public:
    itensor::MPO H0;
    std::array<itensor::MPO,3> M;
    std::array<itensor::MPO,3> M2;
    Hamiltonian H_Details;
    itensor::SiteSet sites;
    itensor::MPO H_flux;

    private:
    double GSE;
    std::vector<std::array<double,2>> E;
    std::vector<std::array<double,2>> Cv;
    std::vector<std::array<double,2>> S;
    std::vector<std::array<double,2>> W;
    std::vector<std::array<double,2>> Chix;
    std::vector<std::array<double,2>> Chiy;
    std::vector<std::array<double,2>> Chiz;
    std::vector<std::array<double,2>> Magx;
    std::vector<std::array<double,2>> Magy;
    std::vector<std::array<double,2>> Magz;


    std::array<int,3> get_neighbour_data_hex(int LX, int LY, int pos);
    std::array<int,3> get_neighbour_data_tri(int LX, int LY, int pos);

    std::array<int,3> get_neighbour_data(int LX, int LY, int pos){
        std::array<int,3> n;
        if (Lattice_Type == 1){
            n = get_neighbour_data_hex(LX,LY,pos);
        } else if (Lattice_Type == 2){
            n = get_neighbour_data_tri(LX,LY,pos);
        }
        return n;
    }


    void add_kitaev_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_magnetic_interaction(int LX, int LY, int aux);
    void add_heisenberg_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gamma_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gammaq_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    itensor::MPO honeycomb_flux_operator(int LX, int LY, int aux);
    std::array<std::array<itensor::MPO,3>,2> magnetization_operators(int LX, int LY, int aux);

    std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M);
    
    void tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, std::array<std::vector<double>,3>& Chi_vec, std::array<std::vector<double>,3>& Mag_vec, itensor::MPS& psi, itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& Sweeps, double& cb);



    template<std::size_t n, typename T>
    void save_data(std::string filename, std::vector<std::array<T,n>>& vec);
    
    template<typename T>
    void save_data(std::string filename, T v);

    std::vector<double> derivative(std::vector<double>& f, double dx);
    std::vector<double> integral(std::vector<double>& f, double dx, double c);
    std::vector<double> multiply(std::vector<double>& a, std::vector<double>& b);



    public:
    template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, T& sites, int aux, std::string shape){
        sites = T((LX*LY+2*aux),{"ConserveQNs=",false});
        auto ampo = itensor::AutoMPO(sites); 
        
        this -> ampo = ampo;
        this -> H_Details = H_Details;
        this -> sites = sites;

        std::vector<int> full_points;
        if (shape == "Honeycomb"){
            Lattice_Type = 1;
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

        
        this -> H0 = itensor::toMPO(this -> ampo);
        auto Ms = magnetization_operators(LX,LY,aux);
        M = Ms[0];
        M2 = Ms[1];

        std::cout << "Spin Type: " << typeid(sites).name() << "\n";
        std::cout << "Lattice Type: " << shape << "\n";
        std::cout << "LX = " << LX << ", LY = " << LY << "\n";
        std::cout << "Auxiliary Sites: " << aux << "\n\n";
        std::cout << "Hamiltonian Parameters \n";
        this->H_Details.print();
        std::cout << "\n";

    }


    void Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, int max_sites=256, int init_rand_sites=32, std::string TDVP_Type="TwoSite");
    std::array<std::vector<std::array<double,2>>,2> Calculate_Heat_Capacity(int TimeSteps, std::vector<double>& intervals, std::vector<std::vector<double>>& Energies);

    Hamiltonian Get_Constants(){
        return H_Details;
    }

    void Print_Interactions(){
        itensor::PrintData(ampo);
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
            std::string xcx = x + "/" + "Chix";
            std::string xcy = x + "/" + "Chiy";
            std::string xcz = x + "/" + "Chiz";
            std::string xmx = x + "/" + "Magx";
            std::string xmy = x + "/" + "Magy";
            std::string xmz = x + "/" + "Magz";
            
            save_data(xE,E);
            save_data(xC,Cv);
            save_data(xS,S);
            save_data(xW,W);
            save_data(xcx,Chix);
            save_data(xcy,Chiy);
            save_data(xcz,Chiz);
            save_data(xmx,Magx);
            save_data(xmy,Magy);
            save_data(xmz,Magz);
        }
        
    }

    
    itensor::MPS DMRG(int Sweeps=10){
        Calc_Type = true;
        auto psi0 = itensor::randomMPS(sites,16);
        auto [energy, psi] = itensor::dmrg(H0,psi0,Sweeps,{"Quiet",true});
        GSE = energy;
        return psi;
    }





};






}










