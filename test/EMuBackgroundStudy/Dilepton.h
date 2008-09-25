#ifndef Dilepton_h
#define Dilepton_h

#include <map>
#include <vector>

class EMuBackgroundStudy;

// DileptonKeys are keys in DileptonMaps to collections of dileptons;
// see the description of DileptonMap below.
class DileptonKey {
 public:
  DileptonKey(unsigned c0, unsigned c1, int q, unsigned ct);

  unsigned collection[2];
  int charge;
  unsigned cuts;
};

// Needed for using DileptonKeys in maps.
class DileptonKeyCmp {
 public:
  bool operator()(const DileptonKey& k1, const DileptonKey& k2);
};

// In this code, a dilepton is just a pair of indices into the tree's
// lepton vectors, and some derived quantities such as mass, delta_r,
// etc.
class Dilepton {
 public:
  Dilepton(EMuBackgroundStudy* events_, unsigned j_, unsigned k_);

  EMuBackgroundStudy* e;
  unsigned j;
  unsigned k;
  double mass;
  double delta_r;
  double delta_phi;
  double d0;
};

// Helper function for sorting by decreasing mass.
bool DileptonDecreasingMass(const Dilepton* a, const Dilepton* b);

// A dilepton collection is a bunch of pointers to dileptons. The
// pointers are owned and managed by EMuBackgroundStudy.
typedef std::vector<Dilepton*> Dileptons;

// Store dileptons accessed by DileptonKeys. The dileptons are grouped
// by which collections of leptons comprise them, what total charge
// they have, and which cuts they passed. See Parameters and
// EMuBackgroundStudy for definitions of these.
typedef std::map<DileptonKey, Dileptons, DileptonKeyCmp> DileptonMap;

#endif // Dilepton_h
