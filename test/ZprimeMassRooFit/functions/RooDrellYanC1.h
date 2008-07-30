/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2007, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROODRELLYANC1
#define ROODRELLYANC1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooDrellYanC1 : public RooAbsPdf {
public:
  RooDrellYanC1() {} ; 
  RooDrellYanC1(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _M,
	      RooAbsReal& _Gamma,
	      RooAbsReal& _theta);
  RooDrellYanC1(const RooDrellYanC1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDrellYanC1(*this,newname); }
  inline virtual ~RooDrellYanC1() { }

protected:

  RooRealProxy x ;
  RooRealProxy M ;
  RooRealProxy Gamma ;
  RooRealProxy theta ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDrellYanC1,1) // Your description goes here...
};
 
#endif
