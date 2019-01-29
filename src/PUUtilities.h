#ifndef Zp2mu_PUUtilities_h
#define Zp2mu_PUUtilities_h

#include <map>
#include <iostream>
#include <TString.h>

namespace PU {

double MC_pileup_weight(int NumTrueInteraction, TString mc, TString data_scale);

}

#endif
