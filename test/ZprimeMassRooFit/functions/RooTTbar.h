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

#ifndef ROOTTBAR
#define ROOTTBAR

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooTTbar : public RooAbsPdf {
public:
  RooTTbar() {} ; 
  RooTTbar(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _alpha,
	      RooAbsReal& _beta,
	      RooAbsReal& _gamma,
	      RooAbsReal& _norm,
	      RooAbsReal& _mpv,
	      RooAbsReal& _sigma);
  RooTTbar(const RooTTbar& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTTbar(*this,newname); }
  inline virtual ~RooTTbar() { }

protected:

  RooRealProxy x ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  RooRealProxy gamma ;
  RooRealProxy norm ;
  RooRealProxy mpv ;
  RooRealProxy sigma ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooTTbar,1) // Your description goes here...
};
 
#endif
