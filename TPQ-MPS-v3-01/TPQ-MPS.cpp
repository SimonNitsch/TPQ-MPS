#include "TPQ-MPS.h"
#pragma once


namespace TPQ_MPS{

std::array<int,3> Kitaev_Model::get_neighbour_data_hex(int LX, int LY, int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (Y != (LY-1)){
        neighbours[0] = pos+1;
    }
    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    if (pos != 1 && pos != LX*LY){
        neighbours[2] = pos-1;
    }
    return neighbours;

}


std::array<int,3> Kitaev_Model::get_neighbour_data_tri(int LX, int LY, int pos){
    int X = (pos-1) / LY;
    int Y = (pos-1) % LY;
    std::array<int,3> neighbours;

    if (X != (LX-1)){
        neighbours[0] = pos + LY;
        if (Y == 0){
            neighbours[2] = pos + 2*LY - 1;
        }
        else{
            neighbours[2] = pos + LY - 1;
        }
    }
    if (Y == (LY-1)){
        neighbours[1] = pos + 1 - LY;
    }
    else{
        neighbours[1] = pos + 1;
    }

    return neighbours;
}





void Kitaev_Model::add_kitaev_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double K = H_Details.get("K");
    
    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += K,"Sx",i+aux,"Sx",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += K,"Sy",i+aux,"Sy",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += K,"Sz",i+aux,"Sz",n[2]+aux;
        }
    }

}


void Kitaev_Model::add_magnetic_interaction(int LX, int LY, int aux){
    double hx = H_Details.get("hx");
    double hy = H_Details.get("hy");
    double hz = H_Details.get("hz");

    for (int i = aux+1; i != LX*LY+aux+1; i++){
        ampo += -1*hx,"Sx",i;
        ampo += -1*hy,"Sy",i;
        ampo += -1*hz,"Sz",i;
    }
}



void Kitaev_Model::add_heisenberg_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double J = H_Details.get("J");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[0]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[0]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[1]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[1]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[2]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[2]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[2]+aux;
        }
    }

}


void Kitaev_Model::add_gamma_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double G = H_Details.get("Gamma");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += G,"Sy",i+aux,"Sz",n[0]+aux;
            ampo += G,"Sz",i+aux,"Sy",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += G,"Sx",i+aux,"Sz",n[1]+aux;
            ampo += G,"Sz",i+aux,"Sx",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += G,"Sx",i+aux,"Sy",n[2]+aux;
            ampo += G,"Sy",i+aux,"Sx",n[2]+aux;
        }
    }

}



void Kitaev_Model::add_gammaq_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double GQ = H_Details.get("GammaQ");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += GQ,"Sy",i+aux,"Sx",n[0]+aux;
            ampo += GQ,"Sx",i+aux,"Sy",n[0]+aux;
            ampo += GQ,"Sz",i+aux,"Sx",n[0]+aux;
            ampo += GQ,"Sx",i+aux,"Sz",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += GQ,"Sx",i+aux,"Sy",n[1]+aux;
            ampo += GQ,"Sy",i+aux,"Sx",n[1]+aux;
            ampo += GQ,"Sz",i+aux,"Sy",n[1]+aux;
            ampo += GQ,"Sy",i+aux,"Sz",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += GQ,"Sx",i+aux,"Sz",n[2]+aux;
            ampo += GQ,"Sz",i+aux,"Sx",n[2]+aux;
            ampo += GQ,"Sz",i+aux,"Sy",n[2]+aux;
            ampo += GQ,"Sy",i+aux,"Sz",n[2]+aux;
        }
    }

}







std::vector<std::array<double,2>> Kitaev_Model::Mean(std::vector<std::vector<double>>& M){
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




void Kitaev_Model::tdvp_loop(std::vector<double>& E_vec, itensor::MPS& psi, itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& Sweeps){
    int count = 0;

    for (int j = 0; j != TimeSteps; j++){
        double E = itensor::tdvp(psi,H0,t,Sweeps,args);
        
        E_vec.push_back(E);
        std::cout << E << "\n";
        
    }
}






template<std::size_t n, typename T>
void Kitaev_Model::save_data(std::string filename, std::vector<std::array<T,n>>& vec){
    std::string ending = ".txt";
    std::string full_file = filename + ending;
    std::ofstream file(full_file);

    for (auto& varr : vec){
        auto v = varr.begin();
        file << *v;
        v++;

        for (; v != varr.end(); v++){
            file << ", " << *v;
        }
        file << "\n";

    }
    file.close();
    std::cout << "Data saved as: " << full_file << "\n" << std::flush;
}


template<typename T>
void Kitaev_Model::save_data(std::string filename, T v){
    std::vector<std::array<T,1>> vec;
    std::array<T,1> v_array = {v};
    vec.push_back(v_array);

    save_data(filename, vec);
}







void Kitaev_Model::Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, std::string Accuracy, bool Heat_Capacity, int init_rand_sites){
    Calc_Type = 2;
    
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);

    std::vector<double> E_vec;
    E_vec.reserve((TimeSteps + 1) * intervals.size());

    itensor::Sweeps Sweeps;
    itensor::Args tdvp_args;
    if (Accuracy == "Fast"){
        Sweeps = itensor::Sweeps(1);
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",1});
    } else if (Accuracy == "Slow") {
        Sweeps = itensor::Sweeps(1,1,1024);
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7});
    } else {
        Sweeps = itensor::Sweeps(1,1,1024);
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",1});
    }


    std::vector<itensor::Cplx> T;
    for (auto & it : intervals){
        itensor::Cplx t = -0.5 * it / static_cast<double>(TimeSteps) * itensor::Cplx_1;
        T.emplace_back(t);
    }
    auto t0 = std::chrono::system_clock::now();
    

    for (int i = 0; i != Evols; i++){
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(sites,init_rand_sites);
        auto t = T.begin();

        std::complex<double> E = itensor::innerC(psi,H0,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));
        tdvp_loop(E_vec,psi,*t,TimeSteps,tdvp_args,Sweeps);
        t++;
        
        for (; t != T.end(); t++){
            E_vec.push_back(0);
            tdvp_loop(E_vec,psi,*t,TimeSteps,tdvp_args,Sweeps);
        }

        
        Energies.push_back(E_vec);
        E_vec.clear();
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << "/" << Evols << ", Time Needed: " << time.count() << " seconds\n" << std::flush;


    }

    if (Heat_Capacity){
        auto CvS = Calculate_Heat_Capacity(TimeSteps,intervals,Energies);
        Cv = CvS[0];
        S = CvS[1];
    }

    E = Mean(Energies);

    auto t3 = std::chrono::system_clock::now();
    auto time_total = std::chrono::duration<double>(t3-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_total);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_total);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_total);
    std::cout << "Finished Imaginary Time Evolution, Time Needed: " << hours.count() << " h" << minutes.count() << " min" << seconds.count() << " sec";
}






std::array<std::vector<std::array<double,2>>,2> Kitaev_Model::Calculate_Heat_Capacity(int TimeSteps, std::vector<double>& intervals, std::vector<std::vector<double>>& Energies){
    std::vector<std::vector<double>> Capacity;
    std::vector<double> C_vec;
    std::vector<std::vector<double>> Entropy;
    std::vector<double> En_vec;
    Calc_Type = 3;

    Capacity.reserve(Energies.size());
    C_vec.reserve(TimeSteps*intervals.size());
    Entropy.reserve(Energies.size());
    En_vec.reserve((TimeSteps-1)*intervals.size());

    
    std::vector<double> d_beta;
    d_beta.reserve(intervals.size());
    for (double& i : intervals){
        d_beta.push_back(i / static_cast<double>(TimeSteps));
    }

    std::vector<std::vector<double>> beta2;
    std::vector<double> one_beta2;
    beta2.reserve(intervals.size());
    one_beta2.reserve(TimeSteps-1);
    std::vector<std::vector<double>> beta_1;
    std::vector<double> one_beta_1;
    beta_1.reserve(intervals.size());
    one_beta_1.reserve(TimeSteps-1);

    for (int i = 0; i != intervals.size(); i++){
        int base_beta = 0;
        for (int k = 0; k != i; k++){
            base_beta += intervals[k];
        }

        for (int j = 0; j != TimeSteps; j++){
            double curr_beta = base_beta+(j)*d_beta[i];
            one_beta2.push_back(curr_beta*curr_beta*(-1));
            if (curr_beta == 0){
                curr_beta = 1;
            }
            one_beta_1.push_back(1/curr_beta);
        }
        beta2.push_back(one_beta2);
        beta_1.push_back(one_beta_1);
        one_beta2.clear();
        one_beta_1.clear();
    }

    for (int i = 0; i != Energies.size(); i++){
        double Int_Const = 0;
        double tobeintlast;

        for (int j = 0; j != intervals.size(); j++){
            auto first = Energies[i].begin() + j*(TimeSteps+1)+1;
            auto last = Energies[i].begin() + (j+1)*(TimeSteps+1);
            auto numbers = std::vector<double>(first,last);
            if (j == 0){
                numbers.insert(numbers.begin(),*(first-1));
            } else {
                numbers.insert(numbers.begin(),*(first-2));
            }

            C_vec.push_back(0);
            En_vec.push_back(0);
            std::vector<double> diff = derivative(numbers,d_beta[j]);
            std::vector<double> c = multiply(diff,beta2[j]);
            C_vec.insert(C_vec.end(),c.begin(),c.end());

            std::vector<double> tobeint = multiply(c,beta_1[j]);
            if (j != 0){
                tobeint.insert(tobeint.begin(),tobeintlast);
            }
            tobeintlast = tobeint.back();
            std::vector<double> en = integral(tobeint,d_beta[j],Int_Const);
            En_vec.insert(En_vec.end(),en.begin(),en.end());
        
            std::cout << "  " << Int_Const << "  ";
            Int_Const += en.back();
            for (auto& nu : tobeint){
                std::cout << nu << " ";
            }
        }
        std::cout << "\n\n";
        Capacity.push_back(C_vec);
        Entropy.push_back(En_vec);
        C_vec.clear();
        En_vec.clear();
    }
    
    

    std::array<std::vector<std::array<double,2>>,2> result = {Mean(Capacity), Mean(Entropy)};
    return result;


}



std::vector<double> Kitaev_Model::derivative(std::vector<double>& f, double dx){
    std::vector<double> df;
    df.reserve(f.size()-1);

    for (auto i = f.begin(); i != (f.end()-1); i++){
        df.push_back((*(i+1) - *i)/dx);
    }
    return df;
}


std::vector<double> Kitaev_Model::integral(std::vector<double>& f, double dx, double c){
    std::vector<double> fdx;
    fdx.reserve(f.size());
    double cur = c;

    for (auto i = f.begin(); i != f.end() - 1; i++){
        cur += *i * dx + (*(i+1) - *i) * dx/2;
        fdx.push_back(cur);
    }
    return fdx;
}



std::vector<double> Kitaev_Model::multiply(std::vector<double>&a, std::vector<double>&b){
    int maxsize = std::max(a.size(),b.size());
    std::vector<double> c;
    c.reserve(maxsize);

    for (int i = 0; i != maxsize; i++){
        c.push_back(a[i]*b[i]);
    }
    return c;
}






}





