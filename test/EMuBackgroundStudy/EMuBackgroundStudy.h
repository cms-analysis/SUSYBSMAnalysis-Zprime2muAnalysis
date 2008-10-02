#ifndef EMuBackgroundStudy_h
#define EMuBackgroundStudy_h

#include <map>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"

#include "Dilepton.h"
#include "Parameters.h"
#include "Utilities.h"

#define MAXSIZE 200

class EMuBackgroundStudy {
 public :
  TTree          *fChain;
  Int_t           fCurrent;

  UChar_t         status;
  UInt_t          evt;
  UInt_t          run;
  Float_t         weight;
  UShort_t        proc_id;
  UChar_t         gen_status;
  Short_t         gen_pdg_id[2];
  Float_t         gen_pt[2];
  Float_t         gen_eta[2];
  Float_t         gen_phi[2];
  Float_t         gen_energy[2];
  Float_t         met_x;
  Float_t         met_y;
  UShort_t        trig_bits;
  Int_t           nleptons;
  UChar_t         coll_status;
  UShort_t        collection[MAXSIZE];
  UShort_t        lep_status[MAXSIZE];
  Short_t         pdg_id[MAXSIZE];
  Float_t         pt[MAXSIZE];
  Float_t         eta[MAXSIZE];
  Float_t         phi[MAXSIZE];
  Float_t         energy[MAXSIZE];
  UShort_t        iso_njets[MAXSIZE];
  UShort_t        iso_ntracks[MAXSIZE];
  Float_t         iso_tk[MAXSIZE];
  Float_t         iso_ecal[MAXSIZE];
  Float_t         iso_hcal[MAXSIZE];
  Float_t         tk_pt[MAXSIZE];
  Float_t         tk_d0[MAXSIZE];
  Float_t         tk_dz[MAXSIZE];
  Float_t         tk_chi2dof[MAXSIZE];
  UShort_t        tk_hits_px[MAXSIZE];
  UShort_t        tk_hits_si[MAXSIZE];
  UShort_t        tk_hits_mu[MAXSIZE];
  UShort_t        tk_hits_px_lost[MAXSIZE];
  UShort_t        tk_hits_si_lost[MAXSIZE];
  UShort_t        tk_hits_mu_lost[MAXSIZE];
  Short_t         lepton_id[MAXSIZE];
  Float_t         h_over_e[MAXSIZE];
  Float_t         delta_phi_in[MAXSIZE];
  Float_t         delta_eta_in[MAXSIZE];
  Float_t         sigma_eta_eta[MAXSIZE];
  Float_t         e_seed_over_p_in[MAXSIZE];
  Int_t           njets;
  Float_t         jet_pt[MAXSIZE];
  Float_t         jet_eta[MAXSIZE];
  Float_t         jet_phi[MAXSIZE];
  Float_t         jet_energy[MAXSIZE];

  TBranch        *b_status;
  TBranch        *b_evt;
  TBranch        *b_run;
  TBranch        *b_weight;
  TBranch        *b_proc_id;
  TBranch        *b_gen_status;
  TBranch        *b_gen_pdg_id;
  TBranch        *b_gen_pt;
  TBranch        *b_gen_eta;
  TBranch        *b_gen_phi;
  TBranch        *b_gen_energy;
  TBranch        *b_met_x;
  TBranch        *b_met_y;
  TBranch        *b_trig_bits;
  TBranch        *b_nleptons;
  TBranch        *b_coll_status;
  TBranch        *b_collection;
  TBranch        *b_lep_status;
  TBranch        *b_pdg_id;
  TBranch        *b_pt;
  TBranch        *b_eta;
  TBranch        *b_phi;
  TBranch        *b_energy;
  TBranch        *b_iso_njets;
  TBranch        *b_iso_ntracks;
  TBranch        *b_iso_tk;
  TBranch        *b_iso_ecal;
  TBranch        *b_iso_hcal;
  TBranch        *b_tk_pt;
  TBranch        *b_tk_d0;
  TBranch        *b_tk_dz;
  TBranch        *b_tk_chi2dof;
  TBranch        *b_tk_hits_px;
  TBranch        *b_tk_hits_si;
  TBranch        *b_tk_hits_mu;
  TBranch        *b_tk_hits_px_lost;
  TBranch        *b_tk_hits_si_lost;
  TBranch        *b_tk_hits_mu_lost;
  TBranch        *b_lepton_id;
  TBranch        *b_h_over_e;
  TBranch        *b_delta_phi_in;
  TBranch        *b_delta_eta_in;
  TBranch        *b_sigma_eta_eta;
  TBranch        *b_e_seed_over_p_in;
  TBranch        *b_njets;
  TBranch        *b_jet_pt;
  TBranch        *b_jet_eta;
  TBranch        *b_jet_phi;
  TBranch        *b_jet_energy;

  EMuBackgroundStudy(std::string fn, int max_ent);
  virtual ~EMuBackgroundStudy();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();

  // JMT
  Long64_t jentry;
  int max_entries;
  std::string out_prefix;
  TFile* out_file;
  unsigned event_cuts;
  Parameters::bkg_id_type bkg_id;
  double weight_final;
  Dileptons all_dileptons;
  DileptonMap dileptons;
  double met;
  int ncleanjets;
  std::vector<std::pair<unsigned,unsigned> > cleanjets;
  bool event_passes[2];

  // Helper functions for the lepton branches.
  bool is_muon(unsigned j);
  bool is_electron(unsigned j);
  int charge(unsigned j);

  void UserInit();
  void Clear();
  void PreDileptonCalcs();
  unsigned CalculateCuts(unsigned j, int k=-1);
  bool SkipEvent();
  void MakeDileptons();
  void PostDileptonCalcs();
  bool InterestingEvent();
  void DumpEvent();
  void FillPlots();
  void WritePlots();
  void DrawPlots();
  void DumpMassIntegrals();
  void Finalize();

  std::string histo_name(std::string name, std::string sub, Parameters::bkg_id_type bkg_name=Parameters::nbkgs);
  std::string histo_name(const DileptonKey& key, std::string sub, Parameters::bkg_id_type bkg_name=Parameters::nbkgs);
  TObject* dynamic_histo(std::string name, std::string sub, const char* title, int nbins, double xlo, double xhi, int ybins=-1, double ylo=-1, double yhi=-1);

  std::vector<Utilities::histo_bkg> GetHistos(const char* base_name);
  void StackMassHistos(const char* base_name, const char* name, const char* tex_abbrev);
  void OverlayHistos(const char* base_name, const char* title, const char* xtitle, const char* name, const int rebin = 0, const double max=-1, const bool move_overflow=false);
  void DrawGenTopo(const char* name, const char* which=0);
};

#endif
