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

 #include "RooDrellYanC2.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 

 ClassImp(RooDrellYanC2) 

 RooDrellYanC2::RooDrellYanC2(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _kappa) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   kappa("kappa","kappa",this,_kappa)
 { 
 } 


 RooDrellYanC2::RooDrellYanC2(const RooDrellYanC2& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   kappa("kappa",this,other.kappa)
 { 
 } 



 Double_t RooDrellYanC2::evaluate() const 
 { 
   if (x > 0)
     return exp(-kappa*x)/x/x;
   else
     return -999;
 }
