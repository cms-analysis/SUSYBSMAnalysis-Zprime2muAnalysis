#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class Zprime2muCompositeCandidatePicker : public edm::EDProducer {
public:
  explicit Zprime2muCompositeCandidatePicker(const edm::ParameterSet&);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // Helper stuff.
  struct reverse_mass_sort {
    bool operator()(const pat::CompositeCandidate& lhs, const pat::CompositeCandidate& rhs) {
      return lhs.mass() > rhs.mass();
    }
  };

  void remove_overlap(pat::CompositeCandidateCollection&) const;
  std::vector<reco::TransientTrack> get_transient_tracks(const pat::CompositeCandidate&) const;

  // Evaluate cuts. Return values are pair<cut decision, variable or
  // variables to embed>. The cut decision should be made whether
  // we're actually going to drop the candidate because of this
  // decision or not (controlled by the cut_on_* variables); that will
  // be handled in the loop in produce.
  std::pair<bool, float>             back_to_back_cos_angle(const pat::CompositeCandidate&) const;
  std::pair<bool, CachingVertex<5> > vertex_constrained_fit(const pat::CompositeCandidate&) const;
  std::pair<bool, float>             dpt_over_pt(const pat::CompositeCandidate&) const;

  // If the variable to embed in the methods above is a simple int or
  // float or is going to be embedded wholesale with the generic
  // userData mechanism, we'll do those explicitly in the loop in
  // produce. Otherwise if it's more complicated, the methods here
  // take care of it.
  void embed_vertex_constrained_fit(pat::CompositeCandidate&, const CachingVertex<5>& vtx) const;

  const edm::InputTag src;
  StringCutObjectSelector<pat::CompositeCandidate> selector;

  const unsigned max_candidates;
  const bool do_remove_overlap;

  const bool cut_on_back_to_back_cos_angle;
  const double back_to_back_cos_angle_min;

  const bool cut_on_vertex_chi2;
  const double vertex_chi2_max;

  const bool cut_on_dpt_over_pt;
  const double dpt_over_pt_max;

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
    vertex_chi2_max(cut_on_vertex_chi2 ? cfg.getParameter<double>("vertex_chi2_max") : 1e99),
    cut_on_dpt_over_pt(cfg.existsAs<double>("dpt_over_pt_max")),
    dpt_over_pt_max(cut_on_dpt_over_pt ? cfg.getParameter<double>("dpt_over_pt_max") : 1e99)
{
  produces<pat::CompositeCandidateCollection>();
}

void Zprime2muCompositeCandidatePicker::remove_overlap(pat::CompositeCandidateCollection& cands) const {
  // For the list of CompositeCandidates, find any that share leptons
  // and remove one of them. The sort order of the input is used to
  // determine which of the pair is to be removed: we keep the first
  // one.

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

std::vector<reco::TransientTrack> Zprime2muCompositeCandidatePicker::get_transient_tracks(const pat::CompositeCandidate& dil) const {
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.

  std::vector<reco::TransientTrack> ttv;
  const size_t n = dil.numberOfDaughters();
  for (size_t i = 0; i < n; ++i) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
    assert(mu);
    const reco::TrackRef& tk = patmuon::getPickedTrack(*mu);
    ttv.push_back(ttkb->build(tk));
  }

  return ttv;
}

std::pair<bool, float> Zprime2muCompositeCandidatePicker::back_to_back_cos_angle(const pat::CompositeCandidate& dil) const {
  // Back-to-back cut to kill cosmics.
  assert(dil.numberOfDaughters() == 2);
  const float cos_angle = dil.daughter(0)->momentum().Dot(dil.daughter(1)->momentum()) / dil.daughter(0)->p() / dil.daughter(1)->p();
  return std::make_pair(cos_angle >= back_to_back_cos_angle_min, cos_angle);
}

std::pair<bool, CachingVertex<5> > Zprime2muCompositeCandidatePicker::vertex_constrained_fit(const pat::CompositeCandidate& dil) const {
  // Loose common vertex chi2 cut.
  assert(dil.numberOfDaughters() == 2);
  if (abs(dil.daughter(0)->pdgId()) != 13 || abs(dil.daughter(1)->pdgId()) != 13)
    return std::make_pair(true, CachingVertex<5>()); // pass objects we don't know how to cut on, i.e. e-mu dileptons

  KalmanVertexFitter kvf(true);
  CachingVertex<5> v = kvf.vertex(get_transient_tracks(dil));

  return std::make_pair(v.isValid() && v.totalChiSquared()/v.degreesOfFreedom() <= vertex_chi2_max, v);
}

void Zprime2muCompositeCandidatePicker::embed_vertex_constrained_fit(pat::CompositeCandidate& dil, const CachingVertex<5>& vtx) const {
  if (!vtx.isValid()) {
    dil.addUserFloat("vertex_chi2", 1e8);
    return;
  }

  dil.addUserFloat("vertex_chi2", vtx.totalChiSquared()/vtx.degreesOfFreedom());

  dil.addUserFloat("vertexX", vtx.position().x());
  dil.addUserFloat("vertexY", vtx.position().y());
  dil.addUserFloat("vertexZ", vtx.position().z());
  dil.addUserFloat("vertexXError", sqrt(vtx.error().cxx()));
  dil.addUserFloat("vertexYError", sqrt(vtx.error().cyy()));
  dil.addUserFloat("vertexZError", sqrt(vtx.error().czz()));

  InvariantMassFromVertex imfv;
  static const double muon_mass = 0.1056583;
  InvariantMassFromVertex::LorentzVector p4 = imfv.p4(vtx, muon_mass);
  Measurement1D mass = imfv.invariantMass(vtx, muon_mass);

  dil.addUserFloat("vertexPX", p4.X());
  dil.addUserFloat("vertexPY", p4.Y());
  dil.addUserFloat("vertexPZ", p4.Z());

  dil.addUserFloat("vertexM",      mass.value());
  dil.addUserFloat("vertexMError", mass.error());
}

std::pair<bool, float> Zprime2muCompositeCandidatePicker::dpt_over_pt(const pat::CompositeCandidate& dil) const {
  // Cut on sigma(pT)/pT to reject grossly mismeasured tracks.
  float dpt_over_pt_largest = -1.;
  const size_t n = dil.numberOfDaughters();
  for (size_t i = 0; i < n; ++i) {
    // Only apply this cut to muons.
    if (abs(dil.daughter(i)->pdgId()) == 13) {
      const reco::CandidateBaseRef& lep = dileptonDaughter(dil, i);
      if (lep.isNonnull()) {
	const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);
	if (mu) {
	  const reco::Track* tk = patmuon::getPickedTrack(*mu).get();
	  if (tk) {
	    const double dpt_over_pt = ptError(tk)/tk->pt();
	    if (dpt_over_pt > dpt_over_pt_largest) dpt_over_pt_largest = dpt_over_pt;
	  }
	}
      }
    }
  }
  return std::make_pair(dpt_over_pt_largest < dpt_over_pt_max, dpt_over_pt_largest);
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
    std::pair<bool, float> cos_angle = back_to_back_cos_angle(*c);
    if (cut_on_back_to_back_cos_angle && !cos_angle.first)
      continue;

    // Loose common vertex chi2 cut.
    std::pair<bool, CachingVertex<5> > vertex = vertex_constrained_fit(*c);
    if (cut_on_vertex_chi2 && !vertex.first)
      continue;

    // Loose cut on sigma(pT)/pT of muon tracks.
    std::pair<bool, float> dpt_over_pt_largest = dpt_over_pt(*c);
    if (cut_on_dpt_over_pt && !dpt_over_pt_largest.first)
      continue;

    // Save the dilepton since it passed the cuts, and store the cut
    // variables and other stuff for use later.
    new_cands->push_back(*c);
    new_cands->back().addUserFloat("cos_angle",   cos_angle.second);
    embed_vertex_constrained_fit(new_cands->back(), vertex.second);
    new_cands->back().addUserFloat("dpt_over_pt", dpt_over_pt_largest.second);
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
