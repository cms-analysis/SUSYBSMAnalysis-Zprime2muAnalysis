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

#ifndef ROODRELLYANC2
#define ROODRELLYANC2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooDrellYanC2 : public RooAbsPdf {
public:
  RooDrellYanC2() {} ; 
  RooDrellYanC2(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _kappa);
  RooDrellYanC2(const RooDrellYanC2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDrellYanC2(*this,newname); }
  inline virtual ~RooDrellYanC2() { }

protected:

  RooRealProxy x ;
  RooRealProxy kappa ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDrellYanC2,1) // Your description goes here...
};
 
#endif
