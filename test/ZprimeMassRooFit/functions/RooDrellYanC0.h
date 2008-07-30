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

#ifndef ROODRELLYANC0
#define ROODRELLYANC0

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooDrellYanC0 : public RooAbsPdf {
public:
  RooDrellYanC0() {} ; 
  RooDrellYanC0(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _M,
	      RooAbsReal& _Gamma,
	      RooAbsReal& _theta);
  RooDrellYanC0(const RooDrellYanC0& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDrellYanC0(*this,newname); }
  inline virtual ~RooDrellYanC0() { }

protected:

  RooRealProxy x ;
  RooRealProxy M ;
  RooRealProxy Gamma ;
  RooRealProxy theta ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDrellYanC0,1) // Your description goes here...
};
 
#endif
