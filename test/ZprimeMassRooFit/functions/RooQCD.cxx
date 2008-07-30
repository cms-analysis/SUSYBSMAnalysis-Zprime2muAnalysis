 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

 // -- CLASS DESCRIPTION [PDF] -- 
 // Your description goes here... 

 #include "Riostream.h" 

 #include "RooQCD.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 

 ClassImp(RooQCD) 

 RooQCD::RooQCD(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _beta1,
                        RooAbsReal& _gamma1,
                        RooAbsReal& _beta2,
                        RooAbsReal& _gamma2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   beta1("beta1","beta1",this,_beta1),
   gamma1("gamma1","gamma1",this,_gamma1),
   beta2("beta2","beta2",this,_beta2),
   gamma2("gamma2","gamma2",this,_gamma2)
 { 
 } 


 RooQCD::RooQCD(const RooQCD& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   beta1("beta1",this,other.beta1),
   gamma1("gamma1",this,other.gamma1),
   beta2("beta2",this,other.beta2),
   gamma2("gamma2",this,other.gamma2)
 { 
 } 



 Double_t RooQCD::evaluate() const
 { 
   Double_t p1 = (x-beta1)*(x-beta1)*exp(-x/gamma1);
   Double_t p2 = exp(-beta2*pow(x,gamma2));

   return p1 * p2; 
 }
