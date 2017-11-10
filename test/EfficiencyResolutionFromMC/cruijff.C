Double_t cruijff( Double_t *x, Double_t * par) {

 // par[0] = normalization
 // par[1] = mean
 // par[2] = sigmaL = sigmaR
 // par[3] = alphaL
 // par[4] = alphaR

 double dx = (x[0]-par[1]) ;
 double sigma = par[2] ; //dx<0 ? par[2]: par[3] ;
 double alpha = dx<0 ? par[3]: par[4] ;
 double f = 2*sigma*sigma + alpha*dx*dx ;
 return par[0] * exp(-dx*dx/f); 

}
