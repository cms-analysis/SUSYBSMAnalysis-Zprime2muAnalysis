#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymmetryHelpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

bool computeFitQuantities(const reco::GenParticleCollection& genParticles,
			  int leptonFlavor, bool internalBremOn,
			  AsymFitData& data) {
  static const bool debug = false;

  HardInteraction hi;
  hi.init(leptonFlavor, true);
  hi.Fill(genParticles);
  if (!hi.IsValid())
    return false;

  // Copy the four-vectors into TLorentzVectors, since our code uses
  // those already
  TLorentzVector v_muq, v_my_dil, v_my_mum, v_my_mup;

  v_muq.SetPxPyPzE(hi.quark->p4().x(), hi.quark->p4().y(), hi.quark->p4().z(), hi.quark->p4().e());
  if (internalBremOn) {
    v_my_mum.SetPxPyPzE(hi.lepMinus->p4().x(), hi.lepMinus->p4().y(), hi.lepMinus->p4().z(), hi.lepMinus->p4().e());
    v_my_mup.SetPxPyPzE(hi.lepPlus ->p4().x(), hi.lepPlus ->p4().y(), hi.lepPlus ->p4().z(), hi.lepPlus ->p4().e());
  } 
  else {
    v_my_mum.SetPxPyPzE(hi.lepMinusNoIB->p4().x(), hi.lepMinusNoIB->p4().y(), hi.lepMinusNoIB->p4().z(), hi.lepMinusNoIB->p4().e());
    v_my_mup.SetPxPyPzE(hi.lepPlusNoIB ->p4().x(), hi.lepPlusNoIB ->p4().y(), hi.lepPlusNoIB ->p4().z(), hi.lepPlusNoIB ->p4().e());
  }

  // The 4-vector for the dimuon
  v_my_dil = v_my_mum + v_my_mup;

  // Store pt, rap, phi, mass, phi_cs, and cos_cs in a structure
  data.cos_true = calcCosThetaTrue(v_muq, v_my_mum, v_my_dil, debug);
  data.pT 	= v_my_dil.Pt();
  data.rapidity = v_my_dil.Rapidity();
  data.phi 	= v_my_dil.Phi();
  data.mass 	= v_my_dil.M(); 
  data.qpL 	= v_muq.Pz();
  data.pL       = v_my_dil.Pz();

  data.cos_cs = calcCosThetaCSAnal(v_my_mum.Pz(), v_my_mum.E(),
				   v_my_mup.Pz(), v_my_mup.E(),
				   v_my_dil.Pt(), v_my_dil.Pz(), 
				   data.mass, debug);
  data.phi_cs = calcPhiCSAnal(v_my_mum.Px(), v_my_mum.Py(),
			      v_my_mup.Px(), v_my_mup.Py(),
			      v_my_dil.Pt(), v_my_dil.Eta(), v_my_dil.Phi(),
			      data.mass, debug);

  // we will later assume the p_z of the dilepton is the p_z of the
  // quark we've mistagged this if the quark was actually going the
  // other direction.  translating from kir_anal.F, we looked at the
  // dilepton formed from adding the muon final-state four-vectors
  // (i.e., we ignored small effect of any brem):
  // data.mistag_true = (v_my_dil.Pz()*hi.quark->p4().z() < 0) ? 1 : 0;
  // but, since we have the true resonance four-vector:
  data.mistag_true = (hi.resonance->p4().z()*hi.quark->p4().z() < 0) ? 1 : 0;

  // alternatively, we've mistagged if the signs of cos_cs and
  // cos_true are different
  data.mistag_cs = (data.cos_cs/data.cos_true > 0) ? 0 : 1;

  // Do some extra calculations to determine if event passed acceptance.
  TLorentzVector v_dil, v_mum, v_mup;
  calc4Vectors(data, v_dil, v_mum, v_mup, debug);
  data.cut_status = diRapAccept(v_dil, v_mum, v_mup);

  return true;
}
