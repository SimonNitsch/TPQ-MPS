#include "TPQ-MPS.h"
#pragma once



using namespace itensor;
namespace TPQ_MPS{


MPS Kitaev_Model::mpo_to_tanmps(MPO& H){
    auto Htan = MPS(sitestan);
    //PrintData(Htan);

    auto H1 = H(1);
    H1 = permute(H1,{rightLinkIndex(H,1),siteInds(H,1)(1),siteInds(H,1)(2)});
    auto ind11 = siteIndex(Htan,1);
    auto indlinkold = sim(rightLinkIndex(H,1));
    auto ten1 = ITensor(ind11,indlinkold);
    //std::cout << "aaaaaaaaaaaaaaaaaaa\n" << std::flush;
    //println(H1);
    //PrintData(ten1);

    for (auto I : iterInds(H1)){
        Cplx value = eltC(H1,I);
        int newind = val(I[1]) + (val(I[2]) - 1) * dims;
        ten1.set(newind,val(I[0]),value);
    }
    Htan.set(1,ten1);
    //std::cout << "aaaaaaaaaaaaaaaaaaa\n" << std::flush;


    for (int i = 2; i != length(Htan); i++){
        auto Hiten = H(i);
        Hiten = permute(Hiten,leftLinkIndex(H,i),rightLinkIndex(H,i),siteInds(H,i)(1),siteInds(H,i)(2));
        auto ind2 = siteInds(Htan,i)(2);
        auto indlinknew = sim(rightLinkIndex(H,i));
        auto ten = ITensor(indlinkold,ind2,indlinknew);
        //PrintData(ten);
        //println(Hiten);

        for (auto I : iterInds(Hiten)){ 
            Cplx value = eltC(Hiten,I);
            int newind = val(I[2]) + (val(I[3]) - 1) * dims;
            ten.set(val(I[0]),newind,val(I[1]),value);
        }
        Htan.set(i,ten);
        indlinkold = indlinknew;
    }

    int last = length(Htan);
    auto Hlast = H(last);
    Hlast = permute(Hlast,{leftLinkIndex(H,last),siteInds(H,last)(1),siteInds(H,last)(2)});

    auto indl2 = siteInds(Htan,last)(2);
    auto tenlast = ITensor(indlinkold,indl2);

    for (auto I : iterInds(Hlast)){
        Cplx value = eltC(Hlast,I);
        int newind = val(I[1]) + (val(I[2]) - 1) * dims;
        tenlast.set(val(I[0]),newind,value);
    }
    Htan.set(last,tenlast);

    return Htan;
}



MPO Kitaev_Model::mpo_to_tanmpo(MPO& H){
    auto Htan = MPO(sitestan);

    auto H1 = H(1);
    H1 = permute(H1,{rightLinkIndex(H,1),siteInds(H,1)(1),siteInds(H,1)(2)});
    auto ind11 = rightLinkIndex(H,1);
    auto ind12 = siteInds(Htan,1)(1);
    auto ind13 = siteInds(Htan,1)(2);
    auto ten1 = ITensor(ind11,ind12,ind13);

    for (auto I : iterInds(H1)){
        Cplx value = eltC(H1,I);        
        for (int j = 0; j != dims; j++){
            int newind1 = val(I[1]) + j * dims;
            int newind2 = val(I[2]) + j * dims;
            ten1.set(val(I[0]),newind1,newind2,value);
        }
    }
    Htan.set(1,ten1);


    for (int i = 2; i != length(Htan); i++){
        auto Hiten = H(i);
        Hiten = permute(Hiten,leftLinkIndex(H,i),rightLinkIndex(H,i),siteInds(H,i)(1),siteInds(H,i)(2));
        auto ind1 = leftLinkIndex(H,i);
        auto ind2 = rightLinkIndex(H,i);
        auto ind3 = siteInds(Htan,i)(1);
        auto ind4 = siteInds(Htan,i)(2);
        auto teni = ITensor(ind1,ind2,ind3,ind4);

        for (auto I : iterInds(Hiten)){
            Cplx value = eltC(Hiten,I);
            for (int j = 0; j != dims; j++){
                int newind1 = val(I[2]) + j * dims;
                int newind2 = val(I[3]) + j * dims;
                teni.set(val(I[0]),val(I[1]),newind1,newind2,value);
            }
        }
        Htan.set(i,teni);

    }

    int last = length(Htan);
    auto Hlast = H(last);
    Hlast = permute(Hlast,{leftLinkIndex(H,last),siteInds(H,last)(1),siteInds(H,last)(2)});

    auto indl1 = leftLinkIndex(H,last);
    auto indl2 = siteInds(Htan,last)(1);
    auto indl3 = siteInds(Htan,last)(2);
    auto tenl = ITensor(indl1,indl2,indl3);

    for (auto I : iterInds(Hlast)){
        Cplx value = eltC(Hlast,I);
        for (int j = 0; j != dims; j++){
            int newind1 = val(I[1]) + j * dims;
            int newind2 = val(I[2]) + j * dims;
            tenl.set(val(I[0]),newind1,newind2,value);
        }
    }
    Htan.set(last,tenl);


    return Htan;   
    
} 


MPS Kitaev_Model::mps_to_tanmps(MPS& X){
    auto Xtan = MPS(sitestan);

    auto X1 = X(1);
    auto indsite1 = siteIndex(Xtan,1);
    auto indold = sim(rightLinkIndex(X,1));
    auto ten1 = ITensor(indsite1,indold);

    for (auto I : iterInds(X1)){
        for (auto J : iterInds(X1)){
            if (val(I[1]) == val(J[1])){   
                auto newind = val(I[0]) + (val(J[0]) - 1) * dims;
                auto value = eltC(X1,I) * eltC(X1,J);
                ten1.set(newind,val(I[1]),value); 
            }
        }
    }
    Xtan.set(1,ten1);

    for (int i = 2; i != length(Xtan); i++){
        auto Xiter = X(i);
        auto indsite = siteInds(Xtan,i)(2);
        auto indlinknew = sim(rightLinkIndex(X,i));
        auto ten = ITensor(indold,indsite,indlinknew);
        for (auto I : iterInds(Xiter)){
            for (auto J : iterInds(Xiter)){
                if (val(I[0]) == val(J[0]) && val(I[2]) == val(I[2])){
                    auto newind = val(I[1]) + (val(J[1]) - 1) * dims;
                    auto value = eltC(Xiter,I) * eltC(Xiter,J);
                    ten.set(val(I[0]),newind,val(I[2]),value);
                }
            }
        }
        Xtan.set(i,ten);
        indold = indlinknew;
    }

    auto last = length(Xtan);
    auto Xlast = X(last);
    auto indsite = siteInds(Xtan,last)(2);
    auto tenlast = ITensor(indold,indsite);
    for (auto I : iterInds(Xlast)){
        for (auto J : iterInds(Xlast)){
            if (val(I[0]) == val(I[0])){
                auto newind = val(I[1]) + (val(J[1]) - 1) * dims;
                auto value = eltC(Xlast,I) * eltC(Xlast,J);
                tenlast.set(val(I[0]),newind,value);
            }
        }
    }
    Xtan.set(last,tenlast);
    return Xtan;


}


double Kitaev_Model::tan_energy(MPS& Hexptan, MPO& Htan, int states, int statesdim){
    double Eprop = 0;
    for (int j = 0; j != states; j++){
        auto ranpsi = randomMPS(sites,statesdim);
        auto ranpsitan = mps_to_tanmps(ranpsi);
        std::complex<double> curE = (innerC(Hexptan,Htan,ranpsitan) / innerC(Hexptan,ranpsitan));
        Eprop += std::real(curE);
    }
    Eprop /= static_cast<double>(states);
    return Eprop;
}



void Kitaev_Model::tan_tdvp_loop(int steps, double dt, MPS& Hexptan, MPO& H0tan, Sweeps& sweeps, Args& tdvp_args, std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, MPO& Hfluxtan, double& cb){

    for (int i = 0; i != steps; i++){
        auto stupid_result = tdvp(Hexptan,H0tan,-1.*dt*Cplx_1,sweeps,tdvp_args);

        cb += dt;
        double current_energy = tan_energy(Hexptan,H0tan);
        double curc = cb * cb * (E_vec.back() - current_energy) / dt;
        double curs = S_vec.back() + dt * 0.5 * (curc/cb + E_vec.back()/(cb-dt));
        double curw = tan_energy(Hexptan,Hfluxtan);

        E_vec.push_back(current_energy);
        C_vec.push_back(curc);
        S_vec.push_back(curs);
        W_vec.push_back(curw);
    }

}






void Kitaev_Model::Tan_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, int max_sites){
    auto Hexp = toExpH(ampo,-1*intervals[0]/static_cast<double>(TimeSteps));
    auto Hexptan = mpo_to_tanmps(Hexp);
    auto H0tan = mpo_to_tanmpo(H0);
    auto Hfluxtan = mpo_to_tanmpo(H_flux);

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = max_sites;
    sweeps.cutoff() = 1e-10;
    sweeps.niter() = 30;

    auto tdvp_args = Args({"Silent",true,"ErrGoal",1E-7,"DoNormalize",false});
    
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);
    std::vector<double> E_vec;
    E_vec.reserve((TimeSteps) * intervals.size() + 1);

    std::vector<std::vector<double>> Capacities;
    Capacities.reserve(Evols);
    std::vector<double> C_vec;
    C_vec.reserve((TimeSteps) * intervals.size() + 1);

    std::vector<std::vector<double>> Entropies;
    Entropies.reserve(Evols);
    std::vector<double> S_vec;
    S_vec.reserve((TimeSteps) * intervals.size() + 1);

    std::vector<std::vector<double>> Flux;
    Flux.reserve(Evols);
    std::vector<double> W_vec;
    W_vec.reserve((TimeSteps) * intervals.size() + 1);



    for (int i = 0; i != Evols; i++){
        E_vec.push_back(0.);
        E_vec.push_back(tan_energy(Hexptan,H0tan));
        W_vec.push_back(0.);
        W_vec.push_back(tan_energy(Hexptan,Hfluxtan));

        double dt1 = intervals[0] / static_cast<double>(TimeSteps);
        double cur_beta = dt1;
        
        C_vec.push_back(0.);
        C_vec.push_back(cur_beta * (E_vec[1] - E_vec[0]));
        S_vec.push_back(0.);
        S_vec.push_back(dt1 * 0.5 * C_vec[1] / cur_beta);


        tan_tdvp_loop(TimeSteps-1,dt1,Hexptan,H0tan,sweeps,tdvp_args,E_vec,C_vec,S_vec,W_vec,Hfluxtan,cur_beta);

        for (auto i = intervals.begin()+1; i != intervals.end(); i++){
            double dt = *i / static_cast<double>(TimeSteps);
            tan_tdvp_loop(TimeSteps,dt,Hexptan,H0tan,sweeps,tdvp_args,E_vec,C_vec,S_vec,W_vec,Hfluxtan,cur_beta);
        }
        Energies.push_back(E_vec);
        E_vec.clear();
        Capacities.push_back(C_vec);
        C_vec.clear();
        Entropies.push_back(S_vec.back() - S_vec);
        S_vec.clear();
        Flux.push_back(W_vec);
        W_vec.clear();

    }

    E = Mean(Energies);
    Cv = Mean(Capacities);
    S = Mean(Entropies);
    W = Mean(Flux);



}






}



