#ifndef CutNrs_h
#define CutNrs_h

namespace cutnrs {
  class HEEPV70 {
  public:
    enum CutIndex {
      ET=0,ETA,DETAINSEED,DPHIIN,SIGMAIETAIETA,E2X5OVER5X5,HADEM,TRKISO,EMHADD1ISO,DXY,MISSHITS,ECALDRIVEN
    };
    
    constexpr static unsigned int kMaxBitNr=ECALDRIVEN;
    constexpr static unsigned int kFullMask=( 0x1 << (kMaxBitNr+1) ) -1;
  };
}

#endif
