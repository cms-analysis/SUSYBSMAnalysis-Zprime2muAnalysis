#ifndef UNBINNEDFITTER
#define UNBINNEDFITTER

//Class to do unbinned maximum likelihood fit.
//Inspired by TTreePlayer::UnbinnedFit in Root's TTreePlayer.cxx
//But this version gets passed arrays of data rather then getting it from tree
//
// Bob Cousins 8/2002

#include <vector>

#include "TVirtualFitter.h"

// Optional variants of extended maximum likelihood for mass reach studies.
// EML_1 lets an overall normalization float; EML_2 introduces two free
// parameters, Nsig and Nbkg, in place of SigFr.
//#define EML_1
//#define EML_2

extern void UnbinnedFitLikelihoodFCN(int &npar, double *gin, double &f,
				     double *u, int flag);

class UnbinnedFitter {
 public:
  UnbinnedFitter();

  TVirtualFitter* getFitter();

  int unbinnedFitExec(const char *funcname,
		      const Option_t *option,	const int nentries,
		      double *data1, double *data2, double *data3,
		      double *weight, double& log_ML);
  int unbinnedFitExec(const char *funcname,
		      const Option_t *option, const int nentries,
		      std::vector<double*> data,
		      double *weight, double& log_ML);
};
#endif
