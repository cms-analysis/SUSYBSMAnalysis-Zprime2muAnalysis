#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/triggerTurnOnElectrons.h"

namespace trigEle_2016{


  TRandom3 randNrGen;
  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5){
    float eff = 0.0;
    eff = 0.5*p0*(1+erf((Et-p1)/(1.414*p2)))+0.5*p3*(1+erf((Et-p4)/(1.414*p5)));
    return eff;
  }
  
  float turnOn(float scEt,float scEta){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.8066, 32.92, 0.5816, 0.1747, 32.66, 1.441);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.901, 33.02, 0.7304, 0.08336, 33.37, 2.464);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.6564, 32.97, 0.7031, 0.3334, 33.33, 1.533);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8954, 33.05, 1.022, 0.0973, 32.32, 4.124);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.7767, 33.08, 0.9646, 0.2196, 33.36, 2.014);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.3626, 34.05, 1.935, 0.6339, 33.59, 1.064);
    else
      return -1.0;
  }
  float turnOn_MW(float scEt,float scEta){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.8315, 33.03, 0.6641, 0.1618, 32.81, 1.666);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.8558, 33.22, 0.757, 0.1391, 33.4, 2.021);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.6471, 33.21, 0.7441, 0.3484, 33.66, 1.642);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.9065, 33.11, 1.079, 0.08481, 31.92, 5.275);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.8312, 33.48, 1.157, 0.1627, 34.04, 2.213);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.9055, 34.45, 1.538, 0.08972, 37.21, 1.672);
    else
      return -1.0;
  }
  
  bool passTrig(float scEt,float scEta){return (11.021/36.459)*turnOn(scEt,scEta)+(25.438/36.459)*turnOn_MW(scEt,scEta)>randNrGen.Uniform(0,1);}

}

namespace trigEle_2017{


  TRandom3 randNrGen;
  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5){
    float eff = 0.0;
    eff = 0.5*p0*(1+erf((Et-p1)/(1.414*p2)))+0.5*p3*(1+erf((Et-p4)/(1.414*p5)));
    return eff;
  }
  float turnOn_MW(float scEt,float scEta){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.873, 33.85, 0.5788, 0.1268, 34.55, 1.852);//rereco all
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.9435, 34.1, 0.792, 0.05638, 34.05, 2.65);//rereco all
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9353, 34.46, 0.9275, 0.06455, 33.7, 2.707);//rereco all
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.5344, 34.44, 0.823, 0.4639, 35.59, 2.119);//rereco all
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.5899, 34.69, 0.8646, 0.409, 35.95, 2.017);//rereco all
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.5655, 34.43, 1.094, 0.432, 36.56, 2.387);//rereco all
    else
      return -1.0;
   } 
  bool passTrig(float scEt,float scEta){return turnOn_MW(scEt,scEta)>randNrGen.Uniform(0,1);}

}


//
//version : New HEEP ID and after using EGamma energy scale and smearing
//Only update for DoubleEle25 Run_all now 
//

namespace trigEle_2018{


  TRandom3 randNrGen;
  float turnOnfunction(float Et, float p0, float p1, float p2, float p3, float p4, float p5){
    float eff = 0.0;
    eff = 0.5*p0*(1+erf((Et-p1)/(1.414*p2)))+0.5*p3*(1+erf((Et-p4)/(1.414*p5)));
    return eff;
  }
  float turnOn_MW(float scEt,float scEta, TString run, bool isDoubleEle25){
    if(isDoubleEle25){
    if(run=="Run_A"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.9233, 25.81, 0.4058, 0.0402, 26.52, 3.212);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.9222, 25.76, 0.476, 0.0433, 26.01, 3.607);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9083, 26, 0.6115, 0.05305, 25.34, 2.279);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.7485, 25.66, 0.6386, 0.2108, 25.64, 2.263);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.8184, 25.92, 0.664, 0.1658, 26.22, 1.921);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.6415, 26.44, 0.7869, 0.3368, 27.64, 1.877);
    else
      return -1.0;
    }
    else if(run=="Run_B"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.9389, 25.98, 0.4151, 0.03921, 26.92, 3.764);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.9232, 25.88, 0.4712, 0.05488, 25.74, 2.616);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9148, 26.01, 0.5986, 0.05494, 25.4, 2.374);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8781, 25.63, 0.7273, 0.08844, 25.96, 3.621);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.7837, 25.86, 0.5784, 0.2012, 26.07, 1.734);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.6676, 26.35, 0.7207, 0.3095, 26.94, 1.781);
    else
      return -1.0;
    }
    else if(run=="Run_C"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.942, 25.74, 0.4395, 0.03665, 26.59, 3.556);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.925, 25.66, 0.5052, 0.05254, 25.37, 3.029);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9091, 25.96, 0.6246, 0.0633, 25.57, 2.443);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.9057, 25.62, 0.7926, 0.06658, 26.5, 5.766);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.9257,25.99,0.7565, 0.06092, 26.62, 2.996);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.3624, 26.5, 0.625, 0.6151, 27.29, 1.411);
    else
      return -1.0;
    }
    else if(run=="Run_D"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.9363, 25.82, 0.4087, 0.04179, 26.21, 3.094);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.9347, 25.8, 0.477, 0.04476, 25.68, 3.444);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9151,26,0.6229,0.05643,25.47,2.507);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8458,25.52,0.6961,0.1248,25.64,3.23);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.8912,25.71,0.6978,0.09515,25.92,2.382);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.9546,26.8,1.055,0.02528,27.21,4.857);
    else
      return -1.0;
    }
    else if(run=="Run_all"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.9359, 25.72, 0.421, 0.03894, 26.34, 3.339);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.9305, 25.72, 0.4978, 0.04524, 25.66, 3.447);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9186, 26.04, 0.6398, 0.0502, 25.39, 2.545);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8628, 25.9, 0.7446, 0.1042, 26.13, 3.608);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.908, 26.04, 0.7248, 0.07583, 26.83, 2.925);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.69, 26.41, 0.8308, 0.2855, 27.28, 1.859 );
    else
      return -1.0;
       }
   else{std::cout<<"wrong run name"<<std::endl;return 1;}
    }

    else{ // for ele33, NOTE: only "Run_all" is correct 
    if(run=="Rereco_B"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.8734, 33.58, 0.4597, 0.08289, 34.62, 2.323);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.8363, 33.68, 0.52, 0.1251, 34.64, 2.275);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.1104, 33.92, 2.565, 0.8681, 33.83, 0.7143);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.759, 34.41, 0.9982, 0.2302, 33.91, 2.368);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.6415, 34.61, 0.7622, 0.3277, 35.18, 2.016);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.7947, 35.0, 1.165, 0.16, 39.38, 1.666);
    else
      return -1.0;
    }
    else if(run=="Rereco_C"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.8983, 33.68, 0.4539, 0.08661, 34.31, 2.081);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.8547, 33.85, 0.5308, 0.134, 34.41, 2.068);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.07665, 33.72, 2.551, 0.9142, 34.14, 0.7229);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.9075, 34.17, 0.8147, 0.08738, 34.82, 3.026);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.8803,34.49,0.8688, 0.116, 34.85, 2.235);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.8673, 34.32, 1.218, 0.1283, 37.06, 1.22);
    else
      return -1.0;
    }
    else if(run=="Rereco_D"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.9139, 33.68, 0.4031, 0.06823, 34.65, 2.454);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,1.056, 33.91, 0.5862, -0.07221, 34.5, 0.01876);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.9141,34.21,0.6963,0.06704,33.64,2.526);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8747,34.21,0.7027,0.103,34.68,2.89);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.9014,34.37,0.6862,0.08579,34.84,2.284);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.8517,34.07,0.8997,0.1356,34.89,2.129);
    else
      return -1.0;
    }
    else if(run=="Run_all"){
    if (0.0<=fabs(scEta) && fabs(scEta)<0.79)
      return turnOnfunction(scEt,0.908, 33.94, 0.4587, 0.0703, 34.71, 2.518);
    else if (0.79<=fabs(scEta) && fabs(scEta)<1.1)
      return turnOnfunction(scEt,0.8776, 33.94, 0.5184, 0.1005, 34.49, 2.409);
    else if (1.1<=fabs(scEta) && fabs(scEta)<1.4442)
      return turnOnfunction(scEt,0.8783, 34.17, 0.686, 0.09349, 34.07, 2.648);
    else if (1.566<=fabs(scEta) && fabs(scEta)<1.70)
      return turnOnfunction(scEt,0.8118, 33.65, 0.7387, 0.158, 33.79, 2.597);
    else if (1.70<=fabs(scEta) && fabs(scEta)<2.1)
      return turnOnfunction(scEt,0.7073, 33.94, 0.6894, 0.2787, 34.12, 1.71);
    else if (2.1<=fabs(scEta) && fabs(scEta)<2.5)
      return turnOnfunction(scEt,0.5525, 34.8, 0.9178, 0.4273, 35.49, 1.888);
    else
      return -1.0;
       }
   else{std::cout<<"wrong run name"<<std::endl;return 1;}
  }//ele33

  }
  
  bool passTrig(float scEt,float scEta, TString run, bool isEle25){return turnOn_MW(scEt,scEta,run, isEle25)>randNrGen.Uniform(0,1);}

}
