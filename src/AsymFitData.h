#ifndef ASYMFITDATA_H
#define ASYMFITDATA_H

enum CUTSTATUS { NOTCUT=0, ETACUT=0x01, PTCUT=0x02 };

struct AsymFitData {
  double pT;
  double pL;
  double qpL;
  double rapidity;
  double mass;
  double phi;
  double cos_cs;
  double phi_cs;
  double cos_true;
  CUTSTATUS cut_status;
  int mistag_true;
  int mistag_cs;
};

std::ostream& operator<< (std::ostream& out, const AsymFitData& data);

#endif
