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

namespace itensor {

class CustomSpinNitschSite;

using CustomSpinNitsch = BasicSiteSet<CustomSpinNitschSite>;

class CustomSpinNitschSite
    {
    Index s;
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
}
