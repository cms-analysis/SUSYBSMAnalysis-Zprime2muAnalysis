
#ifndef DEFINE_DOUBLEELE_Run2
#define DEFINE_DOUBLEELE_Run2

#include "TMath.h"
#include <iostream>
#include "TString.h"
#include "TRandom3.h"

namespace trigEle_2016{


  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5);
  float turnOn(float scEt,float scEta);
  float turnOn_MW(float scEt,float scEta);
  bool passTrig(float scEt,float scEta);

}

namespace trigEle_2017{


  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5);
  float turnOn_MW(float scEt,float scEta);
  bool passTrig(float scEt,float scEta);

}


//
//version : New HEEP ID and after using EGamma energy scale and smearing
//Only update for DoubleEle25 Run_all now 
//

namespace trigEle_2018{


  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5);
  float turnOn_MW(float scEt,float scEta, TString run, bool isDoubleEle25);
  bool passTrig(float scEt,float scEta, TString run, bool isEle25);

}
  

namespace trigEle33l1{


  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5);
  float turnOn_MW(float scEt,float scEta, TString run, bool isDoubleEle25);
  bool passTrig(float scEt,float scEta, TString run, bool isEle25);

}

#endif
