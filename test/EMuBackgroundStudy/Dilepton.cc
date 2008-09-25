#include <algorithm>

#include "EMuBackgroundStudy.h"
#include "Dilepton.h"
#include "Utilities.h"

DileptonKey::DileptonKey(unsigned c0, unsigned c1, int q, unsigned ct)
  : charge(abs(q)),  // Use only absolute charge for now.
    cuts(ct)
{
  collection[0] = c0;
  collection[1] = c1;
}

bool DileptonKeyCmp::operator()(const DileptonKey& k1, const DileptonKey& k2) {
  if      (k1.collection[0] < k2.collection[0]) return true;
  else if (k1.collection[0] > k2.collection[0]) return false;
  else {
    if      (k1.collection[1] < k2.collection[1]) return true;
    else if (k1.collection[1] > k2.collection[1]) return false;
    else {
      if      (k1.charge < k2.charge) return true;
      else if (k1.charge > k2.charge) return false;
      else {
	if      (k1.cuts < k2.cuts) return true;
	else if (k1.cuts > k2.cuts) return false;
	else {
	  return false;
	}
      }
    }
  }
}

Dilepton::Dilepton(EMuBackgroundStudy* events_, unsigned j_, unsigned k_)
  : e(events_), j(j_), k(k_)
{
  // For two leptons with invariant mass near zero,
  // m^2 \approx 2*(E1 E2 - p1 dot p2)
  double m2 = 2*(e->energy[j]*e->energy[k] - e->pt[j]*e->pt[k]*(sinh(e->eta[j])*sinh(e->eta[k]) + cos(e->phi[j] - e->phi[k])));
  // Clamp at zero in case precision limitations force m2 negative.
  mass = m2 < 0 ? 0 : sqrt(m2);

  delta_phi = Utilities::delta_phi(e->phi[j], e->phi[k]);
  delta_r = Utilities::delta_r(e->eta[j], e->eta[k], e->phi[j], e->phi[k]);

  // Store the maximum of the two leptons' track impact parameter..
  d0 = std::max(fabs(e->tk_d0[j]), fabs(e->tk_d0[k]));
}

bool DileptonDecreasingMass(const Dilepton* a, const Dilepton* b) {
  return a->mass > b->mass;
}
