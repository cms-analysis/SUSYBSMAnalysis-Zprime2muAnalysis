#ifndef Zp2mu_PUUtilities1617_h
#define Zp2mu_PUUtilities1617_h

#include <map>
#include <iostream>
#include <TString.h>

namespace PU_2016 {

double MC_pileup_weight(int NumTrueInteraction, int scale, TString dataset);

}
namespace PU_2017 {

double MC_pileup_weight(int NumTrueInteraction, TString mc, TString data_scale);

}

#endif
