#include "TCanvas.h"
#include "TLegend.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"

#include "EMuBackgroundStudy.h"
#include "Parameters.h"
#include "Utilities.h"

using namespace std;

EMuBackgroundStudy::EMuBackgroundStudy(string fn, int max_ent) {
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fn.c_str());
  if (!f) f = new TFile(fn.c_str());
  TTree* tree = (TTree*)gDirectory->Get("Zprime2muBackgrounds/events");
  Init(tree);

  max_entries = max_ent;
  out_prefix = fn.substr(0, fn.find(".root"));
  string::size_type n = out_prefix.find_last_of('/') + 1;
  out_prefix = "out." + out_prefix.substr(n, out_prefix.size() - n + 1);
}

EMuBackgroundStudy::~EMuBackgroundStudy() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t EMuBackgroundStudy::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t EMuBackgroundStudy::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent)
    fCurrent = chain->GetTreeNumber();
  return centry;
}

void EMuBackgroundStudy::Init(TTree *tree) {
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("status", &status, &b_status);
  fChain->SetBranchAddress("evt", &evt, &b_evt);
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
  fChain->SetBranchAddress("proc_id", &proc_id, &b_proc_id);
  fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);
  fChain->SetBranchAddress("gen_pdg_id", gen_pdg_id, &b_gen_pdg_id);
  fChain->SetBranchAddress("gen_pt", gen_pt, &b_gen_pt);
  fChain->SetBranchAddress("gen_eta", gen_eta, &b_gen_eta);
  fChain->SetBranchAddress("gen_phi", gen_phi, &b_gen_phi);
  fChain->SetBranchAddress("gen_energy", gen_energy, &b_gen_energy);
  fChain->SetBranchAddress("met_x", &met_x, &b_met_x);
  fChain->SetBranchAddress("met_y", &met_y, &b_met_y);
  fChain->SetBranchAddress("trig_bits", &trig_bits, &b_trig_bits);
  fChain->SetBranchAddress("nleptons", &nleptons, &b_nleptons);
  fChain->SetBranchAddress("coll_status", &coll_status, &b_coll_status);
  fChain->SetBranchAddress("collection", collection, &b_collection);
  fChain->SetBranchAddress("lep_status", lep_status, &b_lep_status);
  fChain->SetBranchAddress("pdg_id", pdg_id, &b_pdg_id);
  fChain->SetBranchAddress("pt", pt, &b_pt);
  fChain->SetBranchAddress("eta", eta, &b_eta);
  fChain->SetBranchAddress("phi", phi, &b_phi);
  fChain->SetBranchAddress("energy", energy, &b_energy);
  fChain->SetBranchAddress("iso_njets", iso_njets, &b_iso_njets);
  fChain->SetBranchAddress("iso_ntracks", iso_ntracks, &b_iso_ntracks);
  fChain->SetBranchAddress("iso_tk", iso_tk, &b_iso_tk);
  fChain->SetBranchAddress("iso_ecal", iso_ecal, &b_iso_ecal);
  fChain->SetBranchAddress("iso_hcal", iso_hcal, &b_iso_hcal);
  fChain->SetBranchAddress("tk_pt", tk_pt, &b_tk_pt);
  fChain->SetBranchAddress("tk_d0", tk_d0, &b_tk_d0);
  fChain->SetBranchAddress("tk_dz", tk_dz, &b_tk_dz);
  fChain->SetBranchAddress("tk_chi2dof", tk_chi2dof, &b_tk_chi2dof);
  fChain->SetBranchAddress("tk_hits_px", tk_hits_px, &b_tk_hits_px);
  fChain->SetBranchAddress("tk_hits_si", tk_hits_si, &b_tk_hits_si);
  fChain->SetBranchAddress("tk_hits_mu", tk_hits_mu, &b_tk_hits_mu);
  fChain->SetBranchAddress("tk_hits_px_lost", tk_hits_px_lost, &b_tk_hits_px_lost);
  fChain->SetBranchAddress("tk_hits_si_lost", tk_hits_si_lost, &b_tk_hits_si_lost);
  fChain->SetBranchAddress("tk_hits_mu_lost", tk_hits_mu_lost, &b_tk_hits_mu_lost);
  fChain->SetBranchAddress("lepton_id", lepton_id, &b_lepton_id);
  fChain->SetBranchAddress("h_over_e", h_over_e, &b_h_over_e);
  fChain->SetBranchAddress("delta_phi_in", delta_phi_in, &b_delta_phi_in);
  fChain->SetBranchAddress("delta_eta_in", delta_eta_in, &b_delta_eta_in);
  fChain->SetBranchAddress("sigma_eta_eta", sigma_eta_eta, &b_sigma_eta_eta);
  fChain->SetBranchAddress("e_seed_over_p_in", e_seed_over_p_in, &b_e_seed_over_p_in);
  fChain->SetBranchAddress("njets", &njets, &b_njets);
  fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  fChain->SetBranchAddress("jet_energy", jet_energy, &b_jet_energy);
}

bool EMuBackgroundStudy::is_muon(unsigned j) {
  return abs(pdg_id[j]) == 13;
}

bool EMuBackgroundStudy::is_electron(unsigned j) {
  return abs(pdg_id[j]) == 11;
}

int EMuBackgroundStudy::charge(unsigned j) {
  int x = pdg_id[j];
  return -x/abs(x);
}

void EMuBackgroundStudy::UserInit() {
  string fn = out_prefix + ".root";

  out_file = new TFile(fn.c_str(), "recreate");

  Parameters::init();

  gROOT->SetStyle("Plain");
  //gStyle->SetFillColor(0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetMarkerSize(.1);
  gStyle->SetMarkerStyle(8);
  gStyle->SetGridStyle(3);
  gStyle->SetPaperSize(TStyle::kA4);
}

void EMuBackgroundStudy::Clear() {
  event_cuts = 0;
  bkg_id = Parameters::nbkgs;
  weight_final = 0;
  dileptons.clear();
  // Delete the pointers we own.
  for (unsigned i = 0; i < all_dileptons.size(); ++i)
    delete all_dileptons[i];
  all_dileptons.clear();
}

void EMuBackgroundStudy::PreDileptonCalcs() {
  bkg_id = Parameters::bkg_id(proc_id);
  
  // Take the weight from the tree if it is stored, otherwise
  // calculate it, and then apply a K-factor if necessary.
  weight_final = 1;
  if (weight > 0)
    weight_final = weight;
  else
    weight_final = Parameters::partial_weight(bkg_id) * Parameters::int_lumi;
  weight_final *= Parameters::k_factor(bkg_id);

  // Per-event cuts.

  // Require HLT1NonIsoMuon | HLT2NonIsoMuon to have fired.
  if (trig_bits & 3 == 0) event_cuts |= Parameters::trigger;
  
  // Require MET > 20.
  met = sqrt(met_x*met_x + met_y*met_y);
  if (met < 20) event_cuts |= Parameters::met;
  
  // Count the number of clean jets. "Clean" means no electron within
  // delta R of 0.1 (and includes jet pt > 30 and jet abs(eta) < 2.4,
  // done upstream in the ntuple dump).
  ncleanjets = 0;
  cleanjets.clear();
  for (int i = 0; i < njets; ++i) {
    bool jet_is_el = false;
    int j;
    for (j = 0; j < nleptons; ++j) {
      if (collection[j] == Parameters::electrons && Utilities::delta_r(jet_eta[i], eta[j], jet_phi[i], phi[j]) < 0.1) {
	jet_is_el = true;
	break;
      }
    }
    if (!jet_is_el)  ++ncleanjets;

    // Keep track of which jet was cleaned by which electron.
    if (j == nleptons) j = -1;
    cleanjets.push_back(make_pair(i,j));
  }

  if (ncleanjets < 1) event_cuts |= Parameters::njets;
}

unsigned EMuBackgroundStudy::CalculateCuts(unsigned j, int k) {
  // Include the per-event cuts.
  unsigned cuts = event_cuts;

  if (k >= 0) {
    // If k is specified >= 0, we are to check dilepton cuts for this
    // pair of leptons.

    // Delta R cut of 0.1.
    if (Utilities::delta_r(eta[j], eta[k], phi[j], phi[k]) < 0.1)
      cuts |= Parameters::deltar;
  }
  else {
    // Lepton-only cuts for lepton j.

    // pT > 20 and 80 GeV cuts.
    if (pt[j] < 20) cuts |= Parameters::pt;
    if (pt[j] < 80) cuts |= Parameters::pt80;
    
    float s;
    if (is_muon(j)) {
      // pT > 20 and 80 GeV cuts.
      if (pt[j] < 20) cuts |= Parameters::pt20mu;
      if (pt[j] < 80) cuts |= Parameters::pt80mu;
      
      // Tracker sum pT < 10 and 3 GeV isolation cuts.
      if (iso_tk[j] > 10) cuts |= Parameters::isotk;
      if (iso_tk[j] >  3) cuts |= Parameters::isotk3;
      
      // Track impact parameter < 2.5 mm cut.
      if (fabs(tk_d0[j]) > 0.25) cuts |= Parameters::d0;

      // Track number of silicon hits >= 7 cut.
      if (tk_hits_si[j] < 7) cuts |= Parameters::nsihits;
      
      // chi^2/dof < 10 cut.
      if (tk_chi2dof[j] > 10) cuts |= Parameters::chi2dof;

      // Isolation definition according to the top group.
      s = iso_tk[j] + iso_ecal[j] + iso_hcal[j];
    }
    else if (is_electron(j)) {
      // pT > 20 and 80 GeV cuts.
      if (pt[j] < 20) cuts |= Parameters::pt20el;
      if (pt[j] < 80) cuts |= Parameters::pt80el;

      // Track impact parameter < 400 microns cut.
      if (fabs(tk_d0[j]) > 0.04) cuts |= Parameters::d0;
      
      // A muon can brem and the resulting photon can fake an electron
      // (track + energy in the ecal). Cut electrons with delta R to a
      // muon less than 0.1.
      for (int i = 0; i < nleptons; ++i) {
	if (is_muon(i) && collection[i] == collection[j] &&
	    Utilities::delta_r(eta[i], eta[j], phi[i], phi[j]) < 0.1) {
	  cuts |= Parameters::collemu;
	  break;
	}
      }

      // Top group isolation definition.
      s = iso_tk[j];
    }

    // Top group isolation cut.
    if (1/(1 + s/pt[j]) < 0.92) cuts |= Parameters::isos;
  }

  return cuts;
}

bool EMuBackgroundStudy::SkipEvent() {
  return false;
}

void EMuBackgroundStudy::MakeDileptons() {
  // Loop over all unique pairs of leptons.
  for (int j = 0; j < nleptons; ++j) {
    for (int k = j + 1; k < nleptons; ++k) {
      // If the dilepton is not made of the required collections or
      // leptons, skip it.
      //
      // The arrays come pre-ordered, muons before electrons, so in
      // the following when we look for e.g. collections in
      // valid_collections, we don't have to worry about doing the
      // unique ordering.
      int total_charge = abs(charge(j) + charge(k));
      if (Parameters::valid_collection(collection[j], collection[k]) && Parameters::valid_charge(total_charge)) {
	vector<DileptonKey> keys;
	unsigned cut = CalculateCuts(j) | CalculateCuts(k) | CalculateCuts(j, k);
	for (vector<unsigned>::const_iterator cuts = Parameters::valid_cuts.begin(); cuts != Parameters::valid_cuts.end(); ++cuts) {
	  if (!(*cuts & cut))
	    // We will add this dilepton to the set identified by this
	    // key -- i.e. it was made of the required collections and
	    // leptons, and it passed this set of cuts.
	    keys.push_back(DileptonKey(collection[j], collection[k], total_charge, *cuts));
	}
      
	if (keys.size()) {
	  // Allocate just one Dilepton object for this lepton pair.
	  Dilepton* dil = new Dilepton(this, j, k);
	  // Store a pointer to it in the all_dileptons vector, which
	  // we will use to delete the objects at the end of this
	  // event.
	  all_dileptons.push_back(dil);

	  // Also store a pointer to this Dilepton in the collection
	  // of every key for which the dilepton passed the cuts and
	  // was made of the right collections and charge.
	  for (vector<DileptonKey>::const_iterator it = keys.begin(); it != keys.end(); ++it) {
	    if (dileptons.find(*it) == dileptons.end())
	      dileptons[*it] = Dileptons();
	    dileptons[*it].push_back(dil);
	  }
	}
      }
    }
  }
  
  // Sort each collection by decreasing mass, and only keep the
  // highest-mass dilepton.
  for (DileptonMap::iterator it = dileptons.begin(); it != dileptons.end(); ++it) {
    if (it->second.size() > 1) {
      sort(it->second.begin(), it->second.end(), DileptonDecreasingMass);
      it->second.erase(it->second.begin() + 1, it->second.end());
    }
  }
}

void EMuBackgroundStudy::PostDileptonCalcs() {
  //static const DileptonKey pass_key(1, 4, 0, Parameters::pt | Parameters::isotk3 | Parameters::chi2dof );
  static const DileptonKey pass_key(1, 4, 0, Parameters::pt | Parameters::isotk);
  // Event has at least one opposite-sign e-mu dilepton passing cuts.
  event_passes[0] = dileptons.find(pass_key) != dileptons.end();
  // Passes the above and has invariant mass > 200 GeV.
  event_passes[1] = event_passes[0] && dileptons[pass_key].at(0)->mass > 200;
}

bool EMuBackgroundStudy::InterestingEvent() {
  return false;
  return bkg_id == Parameters::ttjets && event_passes[1];
  DileptonKey k(1, 4, 2, Parameters::pt | Parameters::isotk);
  return bkg_id == Parameters::ttjets && dileptons.find(k) != dileptons.end();
  return gen_pdg_id[0]*gen_pdg_id[1] == -143;
  return bkg_id == Parameters::QCD && event_passes[0];
}

void EMuBackgroundStudy::DumpEvent() {
  using namespace Utilities;

  printf("********************************************************************************\n");
  printf("entry: %li/%li run/evt: %u/%u ", long(jentry), long(fChain->GetEntriesFast()), run, evt);
  printf("bkg_id: %s (%i) weight: %f (%f) status: %s\n",
	 Parameters::bkg_name(bkg_id).c_str(), proc_id, weight_final, weight, binary(status).c_str());
  printf("gen_status: %s gen_pdg_ids: (%i %i)\n",
	 binary(gen_status).c_str(), gen_pdg_id[0], gen_pdg_id[1]);
  
  for (unsigned i = 0; i < 2; ++i)
    if (gen_pdg_id[i] != 0)
      printf("  gen_lepton: pdg_id: %3i pt: %8.3f eta: %5.2f phi: %5.2f energy: %8.3f\n",
	     gen_pdg_id[i], gen_pt[i], gen_eta[i], gen_phi[i], gen_energy[i]);

  printf("trig_bits: %s\nmet_x: %8.3f met_y: %8.3f\n", binary(trig_bits).c_str(), met_x, met_y);

  printf("coll_status: %s nleptons: %i\n", binary(coll_status).c_str(), nleptons);
  for (int i = 0; i < nleptons; ++i) {
    printf("#%i: collection: %i lep_status: %s\n  pdg_id: %3i pt: %8.3f eta: %5.2f phi: %5.2f energy: %8.3f\n",
	   i, collection[i], binary(lep_status[i]).c_str(), pdg_id[i], pt[i], eta[i], phi[i], energy[i]);
    printf("  tk_pt: %8.3f tk_d0: %e tk_dz: %5.2f tk_chi2dof: %5.2f\n",
	   tk_pt[i], tk_d0[i], tk_dz[i], tk_chi2dof[i]);
    printf("  tk_hits/lost: px: %i/%i si: %i/%i mu: %i/%i\n",
	   tk_hits_px[i], tk_hits_px_lost[i], tk_hits_si[i], tk_hits_si_lost[i], tk_hits_mu[i], tk_hits_mu_lost[i]);
    printf("  iso_njets: %i iso_ntracks: %i iso_tk: %e iso_ecal: %e iso_hcal: %e\n",
	   iso_njets[i], iso_ntracks[i], iso_tk[i], iso_ecal[i], iso_hcal[i]);
    if (is_electron(i))
      printf("  lepton_id: %i h_over_e: %e delta_phi_in: %e delta_eta_in: %e\n  sigma_eta_eta: %e e_seed_over_p_in: %e\n",
	     lepton_id[i], h_over_e[i], delta_phi_in[i], delta_eta_in[i], sigma_eta_eta[i], e_seed_over_p_in[i]);
    printf("  cut for: %s\n", Parameters::cut_name(CalculateCuts(i)).c_str());

  }

  printf("njets: %i\n", njets);
  for (int i = 0; i < njets; ++i)
    printf("  #%i: pt: %8.3f eta: %5.2f phi: %5.2f energy: %8.3f\n",
	   i, jet_pt[i], jet_eta[i], jet_phi[i], jet_energy[i]);
  printf("ncleanjets: %i (", ncleanjets);
  for (unsigned i = 0; i < cleanjets.size(); ++i)
    printf("%i (%i), ", cleanjets.at(i).first, cleanjets.at(i).second);
  printf(")\n");
  

  unsigned ndil = all_dileptons.size();
  printf("ndileptons: %i nkeys: %i\n", ndil, dileptons.size());
  if (ndil > 0) {
    for (DileptonMap::const_iterator it = dileptons.begin(); it != dileptons.end(); ++it) {
      printf("key: %s   ndileptons: %i\n", Parameters::key_name(it->first).c_str(), it->second.size());
      for (unsigned i = 0; i < it->second.size(); ++i) {
	const Dilepton* dil = it->second.at(i);
	printf("  #%i: leptons: %i,%i mass: %6.2f dR: %6.4f\n", i, dil->j, dil->k, dil->mass, dil->delta_r);
      }
    }
  }
}

string EMuBackgroundStudy::histo_name(string name, string sub, Parameters::bkg_id_type bkg) {
  if (bkg == Parameters::nbkgs) bkg = bkg_id;
  string n = "h_" + Parameters::bkg_name(bkg) + "_" + name;
  if (sub.size()) n += "_" + sub;
  return n;
}

string EMuBackgroundStudy::histo_name(const DileptonKey& key, string sub, Parameters::bkg_id_type bkg) {
  return histo_name(Parameters::key_name(key), sub, bkg);
}

TObject* EMuBackgroundStudy::dynamic_histo(string name, string sub, const char* title, int nbins, double xlo, double xhi, int ybins, double ylo, double yhi) {
  string n = histo_name(name, sub);
  TObject* h = gDirectory->Get(n.c_str());
  if (h) return h;
  if (ybins != -1)
    h = (TObject*)(new TH2F(n.c_str(), title, nbins, xlo, xhi, ybins, ylo, yhi));
  else
    h = (TObject*)(new TH1F(n.c_str(), title, nbins, xlo, xhi));
  ((TH1F*)h)->Sumw2();
  return h;
}

void EMuBackgroundStudy::FillPlots() {
  // In all the dynamic_histo calls, we are allocating a histogram
  // with operator new, but we do not need to keep track of the
  // pointers ourselves. ROOT will allow us access by
  // TDirectory::Get(histo_name), and then it will delete the
  // histograms when the file is closed.

  ((TH1F*)dynamic_histo("proc_id",   "", "process id",           76, 0, 76))->Fill(proc_id);
  ((TH1F*)dynamic_histo("proc_id_w", "", "process id, weighted", 76, 0, 76))->Fill(proc_id, weight_final);

  for (unsigned j = 0; j < 3; j++) {
    // Make three sets of plots: one for all events, one for events
    // with at least one e-mu dilepton passing pT and isolation cuts,
    // and another for events passing the previous and passing M_inv >
    // 200.
    if (j > 0 && !event_passes[j-1])
      continue;

    string which = j == 0 ? "" : (j == 1 ? "ep" : "ep200");

    ((TH1F*)dynamic_histo("met",          which, "missing #slash{E}_{T}",    100, 0, 300))->Fill(met, weight_final);
    ((TH1F*)dynamic_histo("njets",        which, "#jets/event",               10, 0,  10))->Fill(njets, weight_final);
    ((TH1F*)dynamic_histo("ncleanjets",   which, "# cleaned jets/event",      10, 0,  10))->Fill(ncleanjets, weight_final);
    
    double p1 = double(abs(gen_pdg_id[0])), p2 = double(abs(gen_pdg_id[1]));
    if (p1 < p2) {
      double tmp = p1;
      p1 = p2;
      p2 = tmp;
    }
    ((TH2F*)dynamic_histo("gen_topo",     which, "generated event topo",   17, 0, 17, 17, 0, 17))->Fill(p1, p2, weight_final);

    for (DileptonMap::const_iterator it = dileptons.begin(); it != dileptons.end(); ++it) {
      for (unsigned i = 0; i < it->second.size(); ++i) {
	const Dilepton* dil = it->second.at(i);
	((TH1F*)dynamic_histo(Parameters::key_name(it->first), "mass"     + which,   "dilepton mass",             100, 0,    1000))->Fill(dil->mass, weight_final);
	((TH1F*)dynamic_histo(Parameters::key_name(it->first), "deltaR"   + which,   "dilepton dR(l1, l2)",       100, 0,    6.29))->Fill(dil->delta_r, weight_final);
	((TH1F*)dynamic_histo(Parameters::key_name(it->first), "deltaPhi" + which,   "dilepton dPhi(l1, l2)",     100, 0,    3.15))->Fill(dil->delta_phi, weight_final);
	((TH1F*)dynamic_histo(Parameters::key_name(it->first), "d0"       + which,   "lepton tk |d0| (max)",      100, 1e-6,  0.1))->Fill(dil->d0, weight_final);
      }
    }
  }
}

void EMuBackgroundStudy::WritePlots() {
  out_file->Write();
}

vector<Utilities::histo_bkg> EMuBackgroundStudy::GetHistos(const char* base_name) {
  vector<Utilities::histo_bkg> mass_histos;

  vector<Parameters::bkg_id_type>::const_iterator bkg = Parameters::valid_bkg_ids.begin();
  for (unsigned ibkg = 0; bkg != Parameters::valid_bkg_ids.end(); ++bkg, ++ibkg) {
    char buf[1024];
    snprintf(buf, 1024, base_name, Parameters::bkg_name(*bkg).c_str());
    TH1F* h = (TH1F*)gDirectory->Get(buf);
    if (h) {
      string newh = string(h->GetName()) + "_2";
      h = (TH1F*)h->Clone(newh.c_str());
      mass_histos.push_back(make_pair(h, *bkg));
    }
  }

  return mass_histos;
}

void EMuBackgroundStudy::StackMassHistos(const char* base_name, const char* name, const char* tex_abbrev) {
  vector<Utilities::histo_bkg> mass_histos = GetHistos(base_name);
  if (mass_histos.size() == 0) return;
  sort(mass_histos.begin(), mass_histos.end(), Utilities::HistogramIncreasingIntegral);

  char title[1024];
  snprintf(title, 1024, "Opp. sign %s invariant mass (stacked)", tex_abbrev);

  THStack* mass_histos_stacked = new THStack(name, title);

  TLegend* legend = new TLegend(0.77, 0.57, 0.89, 0.89);
  legend->SetFillColor(0);

  for (vector<Utilities::histo_bkg>::const_iterator it = mass_histos.begin(); it != mass_histos.end(); ++it) {
    TH1F* h = it->first;
    Parameters::bkg_id_type b = it->second;

    h->SetFillColor(Parameters::color(b));
    h->SetLineColor(Parameters::color(b));
    h->SetAxisRange(200, 1000);
    h->Rebin();

    mass_histos_stacked->Add(h);
  }

  for (vector<Utilities::histo_bkg>::reverse_iterator it = mass_histos.rbegin(); it != mass_histos.rend(); ++it)
    legend->AddEntry(it->first, Parameters::bkg_name(it->second).c_str());
  
  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  c->SetLogy();

  mass_histos_stacked->SetMaximum(20);
  mass_histos_stacked->Draw("hist");
  legend->Draw("same");

  TAxis* xaxis = mass_histos_stacked->GetXaxis();
  TAxis* yaxis = mass_histos_stacked->GetYaxis();

  double GeV_per_bin = (xaxis->GetXmax() - xaxis->GetXmin())/xaxis->GetNbins();
  xaxis->SetRangeUser(200, 1000);

  char xtitle[128], ytitle[128];
  snprintf(xtitle, 128, "M_{%s} (GeV)", tex_abbrev);
  snprintf(ytitle, 128, "Events/%i GeV/%i pb^{-1}", int(GeV_per_bin), int(Parameters::int_lumi));
  
  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(ytitle);

  char fn[128];
  snprintf(fn, 128, "histos/%s.png", name);
  c->SaveAs(fn);
}

void EMuBackgroundStudy::OverlayHistos(const char* base_name, const char* title, const char* xtitle, const char* name, const int rebin, const double max, const bool move_overflow) {
  vector<Utilities::histo_bkg> mass_histos = GetHistos(base_name);

  TLegend* legend = new TLegend(0.77, 0.57, 0.89, 0.89);
  legend->SetFillColor(0);

  for (vector<Utilities::histo_bkg>::const_iterator it = mass_histos.begin(); it != mass_histos.end(); ++it) {
    TH1F* h = it->first;
    Parameters::bkg_id_type b = it->second;

    legend->AddEntry(h, Parameters::bkg_name(b).c_str());
    
    h->SetStats(0);
    h->SetTitle(title);
    h->SetLineColor(Parameters::color(b));
    h->SetLineWidth(2);

    h->Scale(1./h->Integral());
    h->SetMinimum(0);
    if (max > 0)
      h->SetMaximum(max);

    int r = rebin;
    while (r-- > 0) h->Rebin();

    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle("Arb. units");
    h->GetYaxis()->SetTitleOffset(1.45);

    if (move_overflow) {
      int nbins = h->GetNbinsX();
      h->SetBinContent(nbins, h->GetBinContent(nbins) + h->GetBinContent(nbins+1));
    }
  }

  sort(mass_histos.begin(), mass_histos.end(), Utilities::HistogramDecreasingMaximum);

  TCanvas* c = new TCanvas("c", "c", 1200, 800);

  for (vector<Utilities::histo_bkg>::const_iterator it = mass_histos.begin(); it != mass_histos.end(); ++it) {
    string opt = "hist e";
    if (it != mass_histos.begin())
      opt += " same";
    it->first->Draw(opt.c_str());
  }

  legend->Draw("same");

  char fn[128];
  snprintf(fn, 128, "histos/h_%s.png", name);
  c->SaveAs(fn);  
}

void EMuBackgroundStudy::DrawGenTopo(const char* name, const char* which) {
  TCanvas* c = new TCanvas();
  char buf[128];
  if (which == 0)
    snprintf(buf, 128, "h_%s_gen_topo", name);
  else
    snprintf(buf, 128, "h_%s_gen_topo_%s", name, which);
    
  TH1F* h = (TH1F*)gDirectory->Get(buf);
  if (!h) return;

  h->SetMarkerSize(1.3);
  h->Scale(1/h->Integral());
  h->SetStats(0);
  h->Draw("text");
 
  char fn[128];
  snprintf(fn, 128, "histos/h_gen_topo_%s.png", name);
  c->SaveAs(fn);

  double f_emu   = h->GetBinContent(h->FindBin(13,11));
  double f_etau  = h->GetBinContent(h->FindBin(15,11));
  double f_mutau = h->GetBinContent(h->FindBin(15,13));

  printf("DrawGenTopo sample: %s e/mu: which: %s %.3f e/tau: %.3f mu/tau: %.3f sum: %.3f\n",
	 name, which, f_emu, f_etau, f_mutau, f_emu + f_etau + f_mutau);
}

void EMuBackgroundStudy::DrawPlots() {
  StackMassHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_mass",     "h_elmu_mass_20", "e#mu");
  StackMassHistos("h_%s_bestMuVHEEP_oppSign_isotkpt80_mass",   "h_elmu_mass_80", "e#mu");

  OverlayHistos("h_%s_met",                                          "#slash{E}_{T} by process (overlaid)",        "#slash{E}_{T} (GeV)", "met", 1);
  OverlayHistos("h_%s_njets",                                        "# jets by process (overlaid)",               "# jets",              "njets");
  OverlayHistos("h_%s_ncleanjets",                                   "# cleaned jets by process (overlaid)",       "# jets",              "ncleanjets");
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_deltaR",           "#DeltaR(e,#mu)",                             "#DeltaR",             "deltar", 1, 0.3);
  OverlayHistos("h_%s_bestMuVbestMu_oppSign_isotkpt_deltaR",         "#DeltaR(#mu,#mu)",                           "#DeltaR",             "deltarmumu", 1, 0.3);
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_d0",               "max lepton |d_{0}|",                         "|d_{0}| (cm)",        "d0", 1, 0.3, true);
							          
  OverlayHistos("h_%s_met_ep",                                       "#slash{E}_{T} by process (overlaid)",        "#slash{E}_{T} (GeV)", "met_ep", 1, 0.23);
  OverlayHistos("h_%s_njets_ep",                                     "# jets by process (overlaid)",               "# jets",              "njets_ep");
  OverlayHistos("h_%s_ncleanjets_ep",                                "# cleaned jets by process (overlaid)",       "# jets",              "ncleanjets_ep");
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_deltaRep",         "#DeltaR(e,#mu)",                             "#DeltaR",             "deltar_ep", 1, 0.22);
  OverlayHistos("h_%s_bestMuVbestMu_oppSign_isotkpt_deltaRep",       "#DeltaR(#mu,#mu)",                           "#DeltaR",             "deltarmumu_ep", 1, 0.3);
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_d0ep",             "max lepton |d_{0}|",                         "|d_{0}| (cm)",        "d0_ep", 1, 0.3, true);

  OverlayHistos("h_%s_met_ep200",                                    "#slash{E}_{T} by process (overlaid)",        "#slash{E}_{T} (GeV)", "met_ep200", 1, 0.3);
  OverlayHistos("h_%s_njets_ep200",                                  "# jets by process (overlaid)",               "# jets",              "njets_ep200");
  OverlayHistos("h_%s_ncleanjets_ep200",                             "# cleaned jets by process (overlaid)",       "# jets",              "ncleanjets_ep200");
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_deltaRep200",      "#DeltaR(e,#mu)",                             "#DeltaR",             "deltar_ep200", 1, 0.4);
  OverlayHistos("h_%s_bestMuVbestMu_oppSign_isotkpt_deltaRep200",    "#DeltaR(#mu,#mu)",                           "#DeltaR",             "deltarmumu_ep200", 1, 0.3);
  OverlayHistos("h_%s_bestMuVHEEP_oppSign_isotkpt_d0ep200",          "max lepton |d_{0}|",                         "|d_{0}| (cm)",        "d0_ep200", 1, 0.3, true);

  DrawGenTopo("WW");
  DrawGenTopo("tW");
  DrawGenTopo("ttjets");

  DrawGenTopo("WW", "ep");
  DrawGenTopo("tW", "ep");
  DrawGenTopo("ttjets", "ep");

  DrawGenTopo("WW", "ep200");
  DrawGenTopo("tW", "ep200");
  DrawGenTopo("ttjets", "ep200");
}

void EMuBackgroundStudy::DumpMassIntegrals() {
  string fn = out_prefix + ".integrals";
  FILE* f = fopen(fn.c_str(), "wt"); // I hate ostream..
  fprintf(f, "%20s%10s%30s%10s %7s %7s %7s %7s\n", "","","","",">200", "+/-", ">400", "+/-");

  const double breaks[3] = { 200, 400, 1000 };

  vector<std::pair<unsigned,unsigned> >::const_iterator coll = Parameters::valid_collections.begin();
  for (unsigned icoll = 0; coll != Parameters::valid_collections.end(); ++coll, ++icoll) {
    vector<int>::const_iterator charge = Parameters::valid_charges.begin();
    for (unsigned icharge = 0; charge != Parameters::valid_charges.end(); ++charge, ++icharge) {
      vector<unsigned>::const_iterator cuts = Parameters::valid_cuts.begin();
      for (unsigned icuts = 0; cuts != Parameters::valid_cuts.end(); ++cuts, ++icuts) {
	vector<Parameters::bkg_id_type>::const_iterator bkg = Parameters::valid_bkg_ids.begin();
	for (unsigned ibkg = 0; bkg != Parameters::valid_bkg_ids.end(); ++bkg, ++ibkg) {

	  if (ibkg == 0)
	    fprintf(f, "\n");

	  if (icharge == 0 && icuts == 0 && ibkg == 0)
	    fprintf(f, "%20s", Parameters::collection_name(*coll).c_str());
	  else
	    fprintf(f, "%20s", "");
	  if (icuts == 0 && ibkg == 0)
	    fprintf(f, "%10s", Parameters::charge_name(*charge).c_str());
	  else
	    fprintf(f, "%10s", "");
	  if (ibkg == 0)
	    fprintf(f, "%30s", Parameters::cut_name(*cuts).c_str());
	  else
	    fprintf(f, "%30s", "");
	  fprintf(f, "%10s", Parameters::bkg_name(*bkg).c_str());

	  DileptonKey key(coll->first, coll->second, *charge, *cuts);
	  string name = histo_name(key, "mass", *bkg);
	  TH1F* hist = (TH1F*)gDirectory->Get(name.c_str());
	  if (hist) {
	    // Perform tail integrals starting at each of the breaks above.
	    int binhi = hist->FindBin(1e9);

	    for (unsigned i = 0; i < 2; ++i) {
	      int binlo = hist->FindBin(breaks[i]);
	      double integral = hist->Integral(binlo, binhi);
	      double error2 = 0;
	      
	      for (unsigned j = binlo; j <= unsigned(binhi); ++j) {
		double e = hist->GetBinError(j);
		error2 += e*e;
	      }
	      
	      fprintf(f, " %7.3f %7.3f", integral, sqrt(error2));
	    }
	  }
	  //else
	  //  printf("can't get hist %s!!!\n", name.c_str());


	  fprintf(f, "\n");
	}
      }
    }
  }

  fclose(f);
}

void EMuBackgroundStudy::Finalize() {
  Clear();
  // Closing the file takes care of deleting all those histos we made.
  out_file->Close();
}  
  
void EMuBackgroundStudy::Loop() {
  const bool debug = false;

  UserInit();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (jentry=0; jentry<nentries;jentry++) {
    if (jentry == max_entries) break;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (!debug && jentry % 1000 == 0)
      fprintf(stderr, "\rEntry: %i", jentry);
     
    Clear();
    PreDileptonCalcs();
    if (SkipEvent())
      continue;
    MakeDileptons();
    PostDileptonCalcs();
    if (debug || InterestingEvent())
      DumpEvent();
    FillPlots();
  }

  fprintf(stderr, "\n");

  WritePlots();
  DrawPlots();

  DumpMassIntegrals();
   
  Finalize();
}
