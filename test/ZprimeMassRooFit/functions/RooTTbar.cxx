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

 #include "RooTTbar.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h"

 #include "TMath.h"

 ClassImp(RooTTbar) 

 RooTTbar::RooTTbar(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _alpha,
                        RooAbsReal& _beta,
                        RooAbsReal& _gamma,
                        RooAbsReal& _norm,
                        RooAbsReal& _mpv,
                        RooAbsReal& _sigma) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   gamma("gamma","gamma",this,_gamma),
   norm("norm","norm",this,_norm),
   mpv("mpv","mpv",this,_mpv),
   sigma("sigma","sigma",this,_sigma)
 { 
 } 


 RooTTbar::RooTTbar(const RooTTbar& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   gamma("gamma",this,other.gamma),
   norm("norm",this,other.norm),
   mpv("mpv",this,other.mpv),
   sigma("sigma",this,other.sigma)
 { 
 } 



 Double_t RooTTbar::evaluate() const 
 { 
   Double_t p1 = (x-beta)*(x-beta)*exp(-x/gamma);

   Double_t p2 = TMath::Landau(x,mpv,sigma);

   return alpha*p1 + norm*p2; 
 }
