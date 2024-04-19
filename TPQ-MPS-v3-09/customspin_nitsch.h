//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//

/*

This is a modified version of customspin.h from the ITesnor library
It additionally contains Sx and Sy operators for the CustomSpin class, which are necessary to run the TPQ-MPS scripts

This header file should replace itensor/itensor/mps/sites/customspin.h

*/


#pragma once

#include "itensor/mps/siteset.h"
#include <vector>


itensor::Cplx operator*(std::vector<itensor::Cplx>& a, std::vector<itensor::Cplx>& b){
    itensor::Cplx c = 0;
    for (int i = 0; i != a.size(); i++){
        c += a[i] * b[i];
    }
    return c;
}



std::vector<std::vector<itensor::Cplx>> operator*(std::vector<std::vector<itensor::Cplx>>& a, std::vector<std::vector<itensor::Cplx>>& b){
    std::vector<std::vector<itensor::Cplx>> c;
    std::vector<itensor::Cplx> c_vec;
    std::vector<itensor::Cplx> a_vec;
            std::vector<itensor::Cplx> b_vec;
    int m = a.size();
    int n = a[0].size();
    c.reserve(m);
    c_vec.reserve(n);
    a_vec.reserve(m);
    b_vec.reserve(m);
    

    for (int i = 0; i != m; i++){
        for (int j = 0; j != n; j++){
            a_vec = a[i];
            for (int k = 0; k != m; k++){
                b_vec.push_back(b[j][k]);
            }
            itensor::Cplx cc = a_vec * b_vec;
            c_vec.push_back(cc);
        }
        c.push_back(c_vec);
        a_vec.clear();
        b_vec.clear();
        c_vec.clear();
    }
    return c;
    
}





std::vector<std::vector<itensor::Cplx>> operator+(std::vector<std::vector<itensor::Cplx>>& a, std::vector<std::vector<itensor::Cplx>>& b){
    
    std::vector<std::vector<itensor::Cplx>> c;
    std::vector<itensor::Cplx> c_vec;
    c.reserve(a.size());
    c_vec.reserve(a[0].size());

    for (int i = 0; i != a.size(); i++){
        for (int j = 0; j != a[0].size(); j++){
            c_vec.push_back(a[i][j] + b[i][j]);
        }
        c.push_back(c_vec);
        c_vec.clear();
    }
    return c;
}


std::vector<std::vector<itensor::Cplx>> operator/(std::vector<std::vector<itensor::Cplx>>& a, itensor::Cplx b){
    std::vector<std::vector<itensor::Cplx>> c;
    std::vector<itensor::Cplx> c_vec;
    c.reserve(a.size());
    c_vec.reserve(a[0].size());

    for (auto& i : a){
        for (auto& j : i){
            itensor::Cplx cc = j / b;
            c_vec.push_back(cc);
        }
        c.push_back(c_vec);
        c_vec.clear();
    }
    return c;
}


std::vector<std::vector<itensor::Cplx>> eyemat(int dims){
    std::vector<std::vector<itensor::Cplx>> e;
    std::vector<itensor::Cplx> e_vec;
    e.reserve(dims);
    e_vec.reserve(dims);

    for (int i = 0; i != dims; i++){
        e_vec = std::vector<itensor::Cplx>(dims);
        e_vec[i] = 1;
        e.push_back(e_vec);
    } 

    return e;
}


std::vector<std::vector<itensor::Cplx>> matexp(std::vector<std::vector<itensor::Cplx>>& A, int loops=10){
    int dims = A.size();
    auto res = eyemat(dims);
    auto cur = eyemat(dims);

    for (int i = 1; i != loops+1; i++){
        cur = cur * A;
        cur = cur / i;
        res = res + cur;
    }
    
    return res;
}








namespace itensor {

    



class CustomSpinNitschSite;

using CustomSpinNitsch = BasicSiteSet<CustomSpinNitschSite>;

class CustomSpinNitschSite
    {
    Index s;
    static bool is_called;
    static std::array<std::vector<std::vector<Cplx>>,3> SpinExp;


    
    static std::array<std::vector<std::vector<Cplx>>,3> ExpMatrix(int DoubleSpin){
        int dims = DoubleSpin + 1;
        std::vector<std::vector<Cplx>> Sx;
        std::vector<std::vector<Cplx>> Sy;
        std::vector<std::vector<Cplx>> Sz;
        Sx.reserve(dims);
        Sy.reserve(dims);
        Sz.reserve(dims);
        std::vector<Cplx> svec(dims);
        for (int i = 0; i != dims; i++){
            Sx.push_back(svec);
            Sy.push_back(svec);
            Sz.push_back(svec);
        }

        double DSdouble = static_cast<double>(DoubleSpin);

        for (int i = 0; i != dims-1; i++){
            double sz = static_cast<double>(i) - DSdouble / 2.; 
            double f = std::sqrt(DSdouble/2.*(DSdouble/2.+1.) - sz*(sz+1.));
            
            Cplx X = f/2. * Cplx_1;
            Cplx Y = -f/2. * Cplx_i;

            Sx[i+1][i] = X * M_PI * Cplx_i;
            Sx[i][i+1] = X * M_PI * Cplx_i;
            Sy[i+1][i] = Y * M_PI * Cplx_i;
            Sy[i][i+1] = -Y * M_PI * Cplx_i;
        }
      
        for (int i = 0; i != dims; i++){
            double sz = static_cast<double>(i) - DSdouble / 2.;
            Sz[i][i] = sz * M_PI * Cplx_i;
        }

        auto expx = matexp(Sx);
        auto expy = matexp(Sy);
        auto expz = matexp(Sz);
        for (auto& i : expx){
            for (auto& j : i){
                std::cout << j << " ";
            }
            std::cout << "\n" << std::flush;
        }


        std::array<std::vector<std::vector<Cplx>>,3> result = {expx,expy,expz};
            
        return result;
        }

        

    public:

    CustomSpinNitschSite() { }

    CustomSpinNitschSite(Index I) : s(I) { }
    
    CustomSpinNitschSite(Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);

        auto DSmax = 3;
        if(args.defined("2S")) 
            {
            DSmax = args.getInt("2S");
            }
        else if(args.defined("S")) 
            {
            auto S = args.getReal("S");
            DSmax = std::round(2*S+0.01);
            }
        else
            {
            error("Must pass named args \"2S\" (integer) or \"S\" (real number) to CustomSpinNitsch");
            }

        if(DSmax < 1) error(tinyformat::format("Invalid spin value %d/2 in CustomSpinNitsch",DSmax));
        
        if (is_called == false){
            CustomSpinNitschSite::SpinExp = CustomSpinNitschSite::ExpMatrix(DSmax);
            CustomSpinNitschSite::is_called = true;
        }

        auto tags = TagSet(tinyformat::format("Site,S=%d/2",DSmax));
        if(args.defined("SiteNumber") )
            {
            auto n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        if(conserveQNs)
            {
            if(conserveSz)
                {
                auto qints = Index::qnstorage(DSmax+1);
                for(int n : range(DSmax+1)) 
                    {
                    qints[n] = QNInt(QN({"Sz",static_cast<QNum::qn_t>(2*n-DSmax)}),1);
                    }
                s = Index(std::move(qints),tags);
                }
            else
                {
                s = Index(QN(),DSmax+1,tags);
                }
            }
        else
            {
            if(conserveSz) throw ITError("ConserveSz cannot be true when ConserveQNs=false");
            s = Index(DSmax+1,tags);
            }
        }


    Index
    index() const { return s; }
    
    IndexVal
    state(std::string const& state)
        {
        auto DSmax = dim(s)-1;
        for(auto n : range(DSmax+1))
            {
            if(state == str(2*n-DSmax)) return s(1+n); //state name is "2Sz"
            }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
        }


	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Op = ITensor(dag(s),sP);
        
        auto DSmax = dim(s)-1;

        if(opname == "Sz")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),sz);
                }
            }
        else
        if(opname == "S+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i+1),std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0)));
                }
            }
        else
        if(opname == "S-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i-1),std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0)));
                }
            }

/*

Changes to CustomSpinNitsch by Simon Nitsch
Adding Sx and Sy Operator

*/


	
	else
	if(opname == "Sx")
	    {
	    for (int i=1; i < dim(s); ++i){
		auto sz = -DSmax/2. + i - 1;
		Op.set(s(i),sP(i+1),0.5*std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0)));
		}
	    for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i-1),0.5*std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0)));
                }
	    }
	else
	if(opname == "Sy"){
	    for(int i = 1; i <= dim(s)-1; ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i+1),-0.5*std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0))*Cplx_i);
                }
            for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i-1),0.5*std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0))*Cplx_i);
                }
	    }
    else
    if(opname == "Expx"){
        auto ExMat = SpinExp[0];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                std::complex<double> Ex = ExMat[i][j];
                PrintData(s);

                if (std::abs(Ex) > 1e-6){
                    Op.set(s(j),sP(i),Ex*Cplx_1);
                }
            }
        }
    }
    else
    if(opname == "Expy"){
        auto ExMat= SpinExp[1];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                std::complex<double> Ex = ExMat[i][j];
                PrintData(s);

                if (std::abs(Ex) > 1e-6){
                    Op.set(s(j),sP(i),Ex*Cplx_1);
                }
            }
        }
    }
    else
    if(opname == "Expz"){
        auto ExMat= SpinExp[2];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                std::complex<double> Ex = ExMat[i][j];
                PrintData(s);

                if (std::abs(Ex) > 1e-6){
                    Op.set(s(j),sP(i),Ex*Cplx_1);
                }
            }
        }
    }

	

/*

End of Changes

*/
	
	else
        if(opname == "U+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                Op.set(s(i),sP(i+1),1.0);
                }
            }
        else
        if(opname == "U-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                Op.set(s(i),sP(i-1),1.0);
                }
            }
        else
        if(opname == "Sz2")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),pow(sz, 2.0));
                }
            }
        else
            {
            throw ITError("Operator " + opname + " name not recognized");
            }

        return Op;
        }
    };



std::array<std::vector<std::vector<Cplx>>,3> CustomSpinNitschSite::SpinExp;
bool CustomSpinNitschSite::is_called = false;





}
