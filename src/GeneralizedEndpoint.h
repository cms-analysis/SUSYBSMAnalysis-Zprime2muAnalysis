// * GeneralizedEndpoint.h
#include <iostream>
#include <map>

#ifndef GENERALIZEDENDPOINT_H_
#define GENERALIZEDENDPOINT_H_

class GeneralizedEndpoint {
 public:
   GeneralizedEndpoint();
   virtual ~GeneralizedEndpoint();
   float GeneralizedEndpointPt(float MuonPt, int MuonCharge, float MuonEta, float MuonPhi, int Mode, int verbose=0);
 private:
   std::map<int,std::map<int,float> > _Correction;
   std::map<int,std::map<int,float> > _CorrectionError;   

 };
#endif /* GENERALIZEDENDPOINT_H_ */
