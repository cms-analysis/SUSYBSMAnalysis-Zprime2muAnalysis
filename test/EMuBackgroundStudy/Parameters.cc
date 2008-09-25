#include <algorithm>

#include "Parameters.h"

using namespace std;

const double Parameters::int_lumi = 100; // 100 pb^-1

const unsigned Parameters::muons = 1;     // bestMuons
const unsigned Parameters::electrons = 4; // heepSelector

// Static members must be instantiated.
vector<Parameters::bkg_id_type> Parameters::valid_bkg_ids;
vector<unsigned> Parameters::valid_cuts;
vector<pair<unsigned,unsigned> > Parameters::valid_collections;
vector<int> Parameters::valid_charges;
map<unsigned,string> Parameters::cut_names;
vector<string> Parameters::collection_names;

void Parameters::init() {
  // Oh to have C++0x with STL-container aggregate initializers!

  // Give names to the cut constants.
  cut_names[none]    = "none";
  cut_names[pt]      = "pt";
  cut_names[isotk]   = "isotk";
  cut_names[ptheep]  = "ptheep";
  cut_names[deltar]  = "deltar";
  cut_names[isotk3]  = "isotk3";
  cut_names[trigger] = "trigger";
  cut_names[isos]    = "isos";
  cut_names[met]     = "met";
  cut_names[d0]      = "d0";
  cut_names[nsihits] = "nsihits";
  cut_names[chi2dof] = "chi2dof";
  cut_names[collemu] = "collemu";
  cut_names[njets]   = "njets";
  cut_names[pt20el]  = "pt20el";
  cut_names[pt20mu]  = "pt20mu";
  cut_names[pt80el]  = "pt80el";
  cut_names[pt80mu]  = "pt80mu";

  // Collection names (order is defined by the ntuple dumper).
  collection_names.push_back("GMR");
  collection_names.push_back("bestMu");
  collection_names.push_back("tightEl");
  collection_names.push_back("HEEPVin"); 
  collection_names.push_back("HEEP");

  // Valid backgrounds.
  valid_bkg_ids.push_back(ttjets);
  valid_bkg_ids.push_back(tW);
  valid_bkg_ids.push_back(WW);
  valid_bkg_ids.push_back(Wjets);
  valid_bkg_ids.push_back(WZ);
  valid_bkg_ids.push_back(ZZ);
  valid_bkg_ids.push_back(DY);
  valid_bkg_ids.push_back(QCD);
//valid_bkg_ids.push_back(photonj);
//valid_bkg_ids.push_back(onia);
//valid_bkg_ids.push_back(muX);
//valid_bkg_ids.push_back(elX);
//valid_bkg_ids.push_back(minbias);
//valid_bkg_ids.push_back(DYtt);
//valid_bkg_ids.push_back(Zp);

  // Valid cuts.
//valid_cuts.push_back(none);
//valid_cuts.push_back(pt);
  valid_cuts.push_back(pt | isotk);
  valid_cuts.push_back(ptheep | isotk);
//valid_cuts.push_back(pt | isotk | deltar);
//valid_cuts.push_back(pt | isotk3);
//valid_cuts.push_back(pt20mu | pt80el | isotk);
//valid_cuts.push_back(pt20el | pt80mu | isotk);
//valid_cuts.push_back(pt);
//valid_cuts.push_back(ptheep);
//valid_cuts.push_back(pt | isotk | njets);
  
  // Valid collections.
  valid_collections.push_back(make_pair(muons, muons));
  valid_collections.push_back(make_pair(electrons, electrons));
  valid_collections.push_back(make_pair(muons, electrons));

  // Valid total charges.
  valid_charges.push_back(0);
  valid_charges.push_back(2);
}

bool Parameters::valid_bkg_id(bkg_id_type bkg_id) {
  return find(valid_bkg_ids.begin(), valid_bkg_ids.end(), bkg_id) != valid_bkg_ids.end();
}

bool Parameters::valid_collection(unsigned c0, unsigned c1) {
  return find(valid_collections.begin(), valid_collections.end(), std::make_pair(c0, c1)) != valid_collections.end();
}

bool Parameters::valid_charge(int charge) {
  return find(valid_charges.begin(), valid_charges.end(), charge) != valid_charges.end();
}

bool Parameters::valid_key(const DileptonKey& key) {
  // A key is valid whether the cut is valid or not. (Cuts being valid
  // only affects whether they are included in the dump of the mass
  // integrals.)
  return valid_collection(key.collection[0], key.collection[1]) && valid_charge(key.charge);
}

string Parameters::bkg_name(bkg_id_type bkg) {
  // Since a bkg_id can be only one of the following, don't bother
  // with a map and just use a simple switch to return the name.
  switch (bkg) {
  case Wjets:   return "Wjets";
  case ttjets:  return "ttjets"; 
  case QCD:     return "QCD"; 
  case photonj: return "photonj"; 
  case onia:    return "onia"; 
  case elX:     return "elX"; 
  case muX:     return "muX"; 
  case minbias: return "minbias"; 
  case WW:  	return "WW";
  case WZ:  	return "WZ";
  case ZZ:  	return "ZZ";
  case DY:  	return "DY";
  case Zp:  	return "Zp";
  case tW:  	return "tW";
  case DYtt:    return "DYtt";
  case fantasy: return "fantasy";
  default:      break;
  }
    
  throw error("invalid bkg_id_type in bkg_name");
}

string Parameters::collection_name(pair<unsigned,unsigned> coll) {
  return collection_names[coll.first] + "V" + collection_names[coll.second];
}

string Parameters::charge_name(int charge) {
  if      (charge == 0) return "oppSign";
  else if (charge == 2) return "sameSign";
  else                  return "invalid"; // shouldn't be able to get here...
}

string Parameters::cut_name(unsigned cut, bool shorten) {
  vector<string> names;
  for (unsigned i = 0; i < ncuts; ++i) {
    unsigned c = 1<<i;
    if (cut & c) names.push_back(cut_names[c]);
  }

  sort(names.begin(), names.end());

  string name;
  unsigned nsize = names.size();
  for (unsigned i = 0; i < nsize; ++i) {
    name += names[i];
    if (!shorten && i != nsize - 1) name += ", ";
  }

  return name;
}

string Parameters::key_name(const DileptonKey& key) {
  if (!valid_key(key))
    throw error("tried to get key_name for invalid key!");

  string name;
  name += collection_name(make_pair(key.collection[0], key.collection[1]));
  name += "_";
  name += charge_name(key.charge);
  name += "_";
  name += cut_name(key.cuts, true);
  return name;
}

Parameters::bkg_id_type Parameters::bkg_id(unsigned id) {
  // See https://twiki.cern.ch/twiki/bin/view/CMS/CSA07ProcessId .
  switch (id) {
  case 27: return minbias;
  case 60: return onia;
  case 61: return QCD;
  case 70: return muX;
  case 71: return WW;
  case 72: return WZ;
  case 73: return ZZ;
  case 74: return DY;
  case 75: return Zp;
  case 76: return tW;
  case 77: return DYtt;
  case 58: // fall through
  case 59: return fantasy; // Higgs and Z'->ee. How useful.
  default: break;
  }
    
  if (id >= 0  && id <= 10) return Wjets;
  if (id >= 22 && id <= 26) return ttjets;
  if (id >= 28 && id <= 47) return QCD;
  if (id >= 48 && id <= 57) return photonj;
  if (id >= 62 && id <= 65) return onia;
  if (id >= 66 && id <= 69) return elX;

  throw error("unknown id in Parameters::bkg_id()");
}

int Parameters::color(bkg_id_type bkg) {
  switch (bkg) {
  case ttjets:    return 2;
  case Wjets:     return 3;
  case WW:        return 4;
  case WZ:        return 5;
  case ZZ:        return 6;
  case DY:        return 7;
  case DYtt:      return 29;
  case tW:        return 1;
  case QCD:       return 42;
  case minbias:   return 30;
  case photonj:   return 46;
  case muX:       return 38;
  case elX:       return 49;
  case onia:      return 40;
  case fantasy:   return 0;
  default: break;
  }

  throw error("unknown id in Parameters::color()");
}

double Parameters::partial_weight(bkg_id_type bkg) {
  switch (bkg) {
  case WW:     return 114.3/845261;
  case WZ:     return 49.9/363291;
  case ZZ:     return 16.1/143113;
  case DY:     return 1797./(3152661*478/515.);      // mumu >200: return 2.52/41927*0.587;
  case DYtt:   return 7559./228195; 
  case tW:     return 62./613791;
  default: break;
  }

  throw error("unknown id in Parameters::partial_weight()");
}

double Parameters::k_factor(bkg_id_type bkg) {
  switch (bkg) {
  case ttjets: return 1.87;
  case DY:     return 1.35;
  case DYtt:   return 1.35;
  case Wjets:  return 61/58.2;
  default:     return 1;
  }
}
