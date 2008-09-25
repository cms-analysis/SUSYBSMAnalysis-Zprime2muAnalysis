#ifndef Parameters_h
#define Parameters_h

#include <stdexcept>
#include <string>
#include <utility> // for make_pair
#include <vector>

#include "Dilepton.h"

class Parameters {
 public:
  typedef std::runtime_error error;

  // Codes to take the CSA07 process ids down to a simpler set. Also
  // includes other backgrounds. bkg_id() does the encoding.
  enum bkg_id_type {
    Wjets   = 1<<0,
    ttjets  = 1<<1,
    QCD     = 1<<2,
    photonj = 1<<3,
    onia    = 1<<4,
    elX     = 1<<5,
    muX     = 1<<6,
    minbias = 1<<7,
    WW 	    = 1<<8,
    WZ 	    = 1<<9,
    ZZ 	    = 1<<10,
    DY 	    = 1<<11,
    Zp 	    = 1<<12,
    tW 	    = 1<<13,
    DYtt    = 1<<14,
    fantasy = 1<<15,
    nbkgs   = 16
  };

  // Codes that can be combined bitwise representing cuts implemented.
  enum cut_type {
    none    = 0,
    pt 	    = 1<<0,
    isotk   = 1<<1,
    ptheep  = 1<<2,
    deltar  = 1<<3,
    isotk3  = 1<<4,
    trigger = 1<<5,
    isos    = 1<<6,
    met     = 1<<7,
    d0 	    = 1<<8,
    nsihits = 1<<9,
    chi2dof = 1<<10,
    collemu = 1<<11,
    njets   = 1<<12,
    pt20el  = 1<<13,
    pt20mu  = 1<<14,
    pt80el  = 1<<15,
    pt80mu  = 1<<16,
    ncuts   = 17
  };

  // Since all members are static, we have to initialize them at some
  // point. Dangerous, but I'm too lazy to implement the singleton
  // pattern.
  static void init();

  // Parameters controlling what dileptons are considered.

  static const double int_lumi;    // in pb^-1;
  static const unsigned muons;     // which muon collection to use:
                                   // 0 = globalMuons
                                   // 1 = bestMuons (from cocktail)
  static const unsigned electrons; // which electron collection to use:
                                   // 2 = selectedLayer1Electrons (using E/gamma id)
                                   // 3 = heepSelectorVincent (for comparison with numbers in HEEP note)
                                   // 4 = heepSelector (HEEP electron id)

  static std::vector<bkg_id_type> valid_bkg_ids;                         // Which background ids to include in plots/print-outs.
  static std::vector<unsigned> valid_cuts;                               // Which combinations of cuts to consider.
  static std::vector<std::pair<unsigned,unsigned> > valid_collections;   // Which pairs of collections to consider (e.g.,
                                                                         // dimuons, dielectrons, muon-electron).
  static std::vector<int> valid_charges;                                 // Which total dilepton charges to consider
                                                                         // (0 = opposite-sign, 2 = same-sign).

  // Test functions to determine whether the parameter is valid
  // (according to the collections above).
  static bool valid_bkg_id(bkg_id_type bkg_id);
  static bool valid_collection(unsigned c0, unsigned c1);
  static bool valid_charge(int charge);
  static bool valid_key(const DileptonKey& key);

  // Helpful items for drawing or print-outs.

  // Names for the above cuts.
  static std::map<unsigned, std::string> cut_names;

  // Names for the collections (indexed by number, see "muons" and
  // "electrons" members above).
  static std::vector<std::string> collection_names;

  // Retrieve the names defined above.
  static std::string bkg_name(bkg_id_type bkg);
  static std::string collection_name(std::pair<unsigned,unsigned> coll);
  static std::string charge_name(int charge);
  static std::string cut_name(unsigned cut, bool shorten=false);
  static std::string key_name(const DileptonKey& key);

  // Convert the CSA07 proc id into our own bkg_id, sometimes encoding
  // multiple CSA07 ids into one of ours.
  static bkg_id_type bkg_id(unsigned id);

  // Return a ROOT color number for each defined bkg.
  static int color(bkg_id_type bkg);

  // Calculate the partial weight (sigma * Br / Nevents) if not
  // expecting to find it in the ntuple.
  static double partial_weight(bkg_id_type bkg);

  // K-factors, all to NLO.
  static double k_factor(bkg_id_type bkg);
};

#endif // Parameters_h
