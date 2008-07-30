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

#ifndef ROOQCD
#define ROOQCD

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooQCD : public RooAbsPdf {
public:
  RooQCD() {} ; 
  RooQCD(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _beta1,
	      RooAbsReal& _gamma1,
	      RooAbsReal& _beta2,
	      RooAbsReal& _gamma2);
  RooQCD(const RooQCD& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooQCD(*this,newname); }
  inline virtual ~RooQCD() { }

protected:

  RooRealProxy x ;
  RooRealProxy beta1 ;
  RooRealProxy gamma1 ;
  RooRealProxy beta2 ;
  RooRealProxy gamma2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooQCD,1) // Your description goes here...
};
 
#endif
