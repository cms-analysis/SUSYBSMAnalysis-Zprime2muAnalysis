#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class Zprime2muCompositeCandidatePicker : public edm::EDProducer {
public:
  explicit Zprime2muCompositeCandidatePicker(const edm::ParameterSet&);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  struct reverse_mass_sort {
    bool operator()(const pat::CompositeCandidate& lhs, const pat::CompositeCandidate& rhs) {
      return lhs.mass() > rhs.mass();
    }
  };

  void remove_overlap(pat::CompositeCandidateCollection&) const;
  float back_to_back_cos_angle(const pat::CompositeCandidate&) const;
  std::vector<reco::TransientTrack> get_transient_tracks(const pat::CompositeCandidate&) const;
  float vertex_chi2(const pat::CompositeCandidate&) const;
  void embed_vertex_constrained_fit(pat::CompositeCandidate&) const;

  const edm::InputTag src;
  StringCutObjectSelector<pat::CompositeCandidate> selector;
  const unsigned max_candidates;
  const bool do_remove_overlap;
  const bool cut_on_back_to_back_cos_angle;
  const double back_to_back_cos_angle_min;
  const bool cut_on_vertex_chi2;
  const double vertex_chi2_max;

  edm::ESHandle<TransientTrackBuilder> ttkb;
};

Zprime2muCompositeCandidatePicker::Zprime2muCompositeCandidatePicker(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    selector(cfg.getParameter<std::string>("cut")),
    max_candidates(cfg.getParameter<unsigned>("max_candidates")),
    do_remove_overlap(cfg.getParameter<bool>("do_remove_overlap")),
    cut_on_back_to_back_cos_angle(cfg.existsAs<double>("back_to_back_cos_angle_min")),
    back_to_back_cos_angle_min(cut_on_back_to_back_cos_angle ? cfg.getParameter<double>("back_to_back_cos_angle_min") : -2),
    cut_on_vertex_chi2(cfg.existsAs<double>("vertex_chi2_max")),
    vertex_chi2_max(cut_on_vertex_chi2 ? cfg.getParameter<double>("vertex_chi2_max") : 1e99)
{
  produces<pat::CompositeCandidateCollection>();
}

void Zprime2muCompositeCandidatePicker::remove_overlap(pat::CompositeCandidateCollection& cands) const {
  // Don't bother doing anything if there's just one candidate.
  if (cands.size() < 2) return; 

  pat::CompositeCandidateCollection::iterator p, q;
  for (p = cands.begin(); p != cands.end() - 1; ) {
    for (q = p + 1; q != cands.end(); ++q) {         
      // Check to see if any of the leptons in p is in q also. If so,
      // remove q (e.g. the one with lower invariant mass since we
      // have sorted the vector already), reset pointers and restart.

      // To do this we need the unique ids of the daughters, i.e. the
      // refs into the original lepton collections.
      typedef std::vector<reco::CandidateBaseRef> refs;
      refs prefs, qrefs;
      for (size_t i = 0; i < p->numberOfDaughters(); ++i)
	prefs.push_back(p->daughter(i)->masterClone());
      for (size_t i = 0; i < q->numberOfDaughters(); ++i)
	qrefs.push_back(q->daughter(i)->masterClone());

      // Compare every pair of (pref, qref) to check for any lepton
      // being shared.
      bool any_shared = false;
      refs::const_iterator pr = prefs.begin(), pre = prefs.end(),
	qr = qrefs.begin(), qre = qrefs.end();
      for ( ; pr != pre && !any_shared; ++pr)
	for ( ; qr != qre && !any_shared; ++qr)
	  if (pr == qr)
	    any_shared = true;

      if (any_shared) {
	cands.erase(q);
	p = cands.begin();
      }
      else
	++p;
    }
  }
}

float Zprime2muCompositeCandidatePicker::back_to_back_cos_angle(const pat::CompositeCandidate& dil) const {
  assert(dil.numberOfDaughters() == 2);
  return dil.daughter(0)->momentum().Dot(dil.daughter(1)->momentum()) / dil.daughter(0)->p() / dil.daughter(1)->p();
}

std::vector<reco::TransientTrack> Zprime2muCompositeCandidatePicker::get_transient_tracks(const pat::CompositeCandidate& dil) const {
  std::vector<reco::TransientTrack> ttv;

  const size_t n = dil.numberOfDaughters();
  for (size_t i = 0; i < n; ++i) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
    if (mu == 0)
      continue;
    const reco::TrackRef& tk = patmuon::getPickedTrack(*mu);
    ttv.push_back(ttkb->build(tk));
  }

  return ttv;
}

float Zprime2muCompositeCandidatePicker::vertex_chi2(const pat::CompositeCandidate& dil) const {
  assert(dil.numberOfDaughters() == 2);
  if (abs(dil.daughter(0)->pdgId()) != 13 || abs(dil.daughter(1)->pdgId()) != 13)
    return -999; // pass objects we don't know how to cut on, i.e. e-mu dileptons

  KalmanVertexFitter kvf(true);
  TransientVertex tv = kvf.vertex(get_transient_tracks(dil));

  return tv.isValid() ? tv.normalisedChiSquared() : 1e8;
}

void Zprime2muCompositeCandidatePicker::embed_vertex_constrained_fit(pat::CompositeCandidate& dil) const {
  assert(dil.numberOfDaughters() == 2);
  if (abs(dil.daughter(0)->pdgId()) != 13 || abs(dil.daughter(1)->pdgId()) != 13)
    return; // don't do anything if it's not a dimuon

  KinematicParticleFactoryFromTransientTrack pfactory;
  std::vector<RefCountedKinematicParticle> muons;
  std::vector<reco::TransientTrack> ttv = get_transient_tracks(dil);
  const size_t n = ttv.size();
  assert(n >= 2);
  float m_sigma = 0.0000001;
  for (size_t i = 0; i < n; ++i)
    muons.push_back(pfactory.particle(ttv[i], 0.1056583, 0., 0., m_sigma));

  KinematicParticleVertexFitter fitter;
  RefCountedKinematicTree tree = fitter.fit(muons);

  if (tree->isValid()) {
    tree->movePointerToTheTop();
    RefCountedKinematicParticle fitted_dimu = tree->currentParticle();
    const AlgebraicVector7& params = fitted_dimu->currentState().kinematicParameters().vector();
    const AlgebraicSymMatrix77& covmat = fitted_dimu->currentState().kinematicParametersError().matrix();

    static const char* names[7] = { "X", "Y", "Z", "PX", "PY", "PZ", "M" };
    char buf[256];
    for (int i = 0; i < 7; ++i) {
      snprintf(buf, 256, "vertex%s", names[i]);
      dil.addUserFloat(buf, params(i));
      snprintf(buf, 256, "vertex%sError", names[i]);
      dil.addUserFloat(buf, sqrt(covmat(i,i)));
    }
  }
}

void Zprime2muCompositeCandidatePicker::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<pat::CompositeCandidateCollection> cands;
  event.getByLabel(src, cands);
  
  // does this get cached correctly? do we care?
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);

  std::auto_ptr<pat::CompositeCandidateCollection> new_cands(new pat::CompositeCandidateCollection);

  // Copy all the candidates that pass the specified cuts into the new
  // output vector. Also embed into the output dimuons any other
  // things that are best to just calculate once and for all.
  for (pat::CompositeCandidateCollection::const_iterator c = cands->begin(), ce = cands->end(); c != ce; ++c) {
    // Some cuts can be simply specified through the
    // StringCutSelector.
    if (!selector(*c))
      continue;
   
    // Now apply cuts that can't.

    // Back-to-back cut to kill cosmics.
    const float cos_angle = back_to_back_cos_angle(*c);
    if (cut_on_back_to_back_cos_angle && cos_angle < back_to_back_cos_angle_min)
      continue;

    // Loose common vertex chi2 cut.
    const float vtx_chi2 = vertex_chi2(*c);
    if (cut_on_vertex_chi2 && vtx_chi2 > vertex_chi2_max)
      continue;

    // Save the dilepton since it passed the cuts, and store the cut
    // variables for use later.
    new_cands->push_back(*c);
    new_cands->back().addUserFloat("cos_angle",   cos_angle);
    new_cands->back().addUserFloat("vertex_chi2", vtx_chi2);

    // Add any other userData.
    embed_vertex_constrained_fit(new_cands->back());
  }

  // Sort candidates so we keep the ones with larger invariant
  // mass. Could make configurable to choose other sorting.
  sort(new_cands->begin(), new_cands->end(), reverse_mass_sort());

  // Remove cands of lower invariant mass that are comprised of a
  // lepton that has been used by a higher invariant mass one.
  if (do_remove_overlap)
    remove_overlap(*new_cands);

  // Only return the maximum number of candidates specified.
  if (new_cands->size() > max_candidates)
    new_cands->erase(new_cands->begin() + max_candidates, new_cands->end());
  
  event.put(new_cands);
}

DEFINE_FWK_MODULE(Zprime2muCompositeCandidatePicker);
