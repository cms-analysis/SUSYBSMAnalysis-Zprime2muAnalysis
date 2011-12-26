/** \class Zprime2muAnalysisPlots
 *  Make inclusive kinematic and quality plots
 *
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

// for selection cut
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

#include <TH1.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TDirectory.h>

#include <map>
#include <string>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

class Zprime2muAnalysisPlots: public edm::EDAnalyzer {
    public:
        /// Constructor
        Zprime2muAnalysisPlots(const edm::ParameterSet& pset) ;

        /// Destructor
        virtual ~Zprime2muAnalysisPlots() ;

        // Operations
        void analyze(const edm::Event & event, const edm::EventSetup& eventSetup) ;

        void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);

        void book(const TFileDirectory &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void book(const TFileDirectory &fs, const edm::ParameterSet &pset, const std::string &name) { book(fs,pset,name,name); }

        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { bookProf(fs,pset,name,name); }

    private:
  //edm::InputTag muons_;
        edm::InputTag dilepton_src_;
        StringCutObjectSelector<pat::Muon> selector_;

        edm::InputTag primaryVertices_;
        edm::InputTag normalization_;

        // we don't care too much about performance
        std::map<std::string, std::map<std::string, TH1*> >  plots;
        std::map<std::string, TProfile*> profiles;

        TH1D *luminosity;
};

/// Constructor
Zprime2muAnalysisPlots::Zprime2muAnalysisPlots(const edm::ParameterSet& pset):
    dilepton_src_(pset.getParameter<edm::InputTag>("dilepton_src")),
    selector_(pset.getParameter<std::string>("selection")),
    primaryVertices_(pset.getParameter<edm::InputTag>("primaryVertices")),
    luminosity(0) // by default, we don't have luminosity info
{
    edm::Service<TFileService> fs;

    TFileDirectory md = fs->mkdir("metadata");
    md.cd(); // JMTBAD should check return value...
    TDirectory *md_dir = md.getBareDirectory();
    md_dir->WriteTObject(new TObjString(dilepton_src_.encode().c_str()), "muons");
    md_dir->WriteTObject(new TObjString(pset.getParameter<std::string>("selection").c_str()), "selection");
    
    //
    // make the directories for the daughters
    //
    std::vector<TFileDirectory> dirs;

    TFileDirectory dauDir1 = fs->mkdir("muon1");
    dirs.push_back(dauDir1);
    TFileDirectory dauDir2 = fs->mkdir("muon2");
    dirs.push_back(dauDir2);

    for(std::vector<TFileDirectory>::const_iterator iDir = dirs.begin(); iDir != dirs.end(); ++iDir) {
    book(*iDir, pset, "p"); 
    book(*iDir, pset, "pt"); 
    book(*iDir, pset, "eta"); 
    book(*iDir, pset, "phi"); 
    book(*iDir, pset, "charge"); 

    book(*iDir, pset, "pSta",   "p"); 
    book(*iDir, pset, "ptSta",  "pt"); 
    book(*iDir, pset, "etaSta", "eta"); 
    book(*iDir, pset, "phiSta", "phi"); 

    book(*iDir, pset, "dxyCoarse");
    book(*iDir, pset, "dxyFine");
    book(*iDir, pset, "dzCoarse");
    book(*iDir, pset, "dzFine");

    book(*iDir, pset, "pixelHits");
    book(*iDir, pset, "pixelLayers");
    book(*iDir, pset, "trackerHits");
    book(*iDir, pset, "trackerLostHitsInner",  "trackerLostHits");
    book(*iDir, pset, "trackerLostHitsMiddle", "trackerLostHits");
    book(*iDir, pset, "trackerLostHitsOuter",  "trackerLostHits");
    book(*iDir, pset, "muonHits");
    book(*iDir, pset, "muonBadHits");
    book(*iDir, pset, "globalHits");
    book(*iDir, pset, "globalMuonHits","muonHits");
    book(*iDir, pset, "trackerChi2n");
    book(*iDir, pset, "muonChi2n");
    book(*iDir, pset, "globalChi2n");

    book(*iDir, pset, "trackIso05", "isolation");
    book(*iDir, pset, "ecalIso05",  "isolation");
    book(*iDir, pset, "hcalIso05",  "isolation");
    book(*iDir, pset, "trackIso03", "isolation");
    book(*iDir, pset, "ecalIso03",  "isolation");
    book(*iDir, pset, "hcalIso03",  "isolation");
    book(*iDir, pset, "combRelIso03", "relIso");
    book(*iDir, pset, "combRelIso05", "relIso");

    book(*iDir, pset, "muonStationsValid",    "muonStations");
    book(*iDir, pset, "muonStationsAny",      "muonStations");
    book(*iDir, pset, "muonStationsDTValid",  "muonStations");
    book(*iDir, pset, "muonStationsDTAny",    "muonStations");
    book(*iDir, pset, "muonStationsCSCValid", "muonStations");
    book(*iDir, pset, "muonStationsCSCAny",   "muonStations");
    book(*iDir, pset, "muonStationsRPCValid", "muonStations");
    book(*iDir, pset, "muonStationsRPCAny",   "muonStations");
    book(*iDir, pset, "segmentMatchesArb",     "segmentMatches"); 
    book(*iDir, pset, "segmentMatchesNoArb",   "segmentMatches"); 
    book(*iDir, pset, "segmentMatchesFailArb", "segmentMatches"); 
    book(*iDir, pset, "segmentCompatArb",      "segmentCompat"); 
    book(*iDir, pset, "segmentCompatNoArb",    "segmentCompat"); 
    book(*iDir, pset, "caloCompat",            "caloCompat"); 
    }

    //
    // make the directory for the diLepton
    //
    TFileDirectory diLepDir = fs->mkdir("diLepton");
    std::vector<TFileDirectory> diDirs;
    diDirs.push_back(diLepDir);
    for(std::vector<TFileDirectory>::const_iterator iDir = diDirs.begin(); iDir != diDirs.end(); ++iDir) {
      book(*iDir, pset, "p"); 
      book(*iDir, pset, "pt"); 
      book(*iDir, pset, "eta"); 
      book(*iDir, pset, "phi"); 
      book(*iDir, pset, "charge");
      book(*iDir, pset, "mass");  
    }
    
    if (pset.existsAs<edm::InputTag>("normalization")) {
        normalization_ = pset.getParameter<edm::InputTag>("normalization");
        luminosity = fs->make<TH1D>("normalization", "normalization", 1, 0, 1);
        luminosity->Sumw2();
    }

}

/// Destructor
Zprime2muAnalysisPlots::~Zprime2muAnalysisPlots()
{

}

void Zprime2muAnalysisPlots::book(const TFileDirectory &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    fs.cd(); // JMTBAD should check return value
    TDirectory *fd = fs.getBareDirectory();
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        plots[fd->GetName()][name] = fs.make<TH1D>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        plots[fd->GetName()][name] = fs.make<TH1D>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}

void Zprime2muAnalysisPlots::bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}


void Zprime2muAnalysisPlots::analyze(const edm::Event & event, const edm::EventSetup& eventSetup){
    using namespace edm;
    using namespace std;
    
    Handle<vector<reco::Vertex> > vertices;
    event.getByLabel(primaryVertices_, vertices);

    edm::Handle<reco::CompositeCandidateView> dileptons;
    event.getByLabel(dilepton_src_, dileptons);
    
    if (!dileptons.isValid())
      edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src_ << " and failed!";
    else {
      //edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dilepton_src_ << " and succeeded! " << dileptons->size();
    }

    foreach (const pat::CompositeCandidate & patCompCand, *dileptons) {
      // first fill the diLepton histograms
      {
	// ultimately, I would like to use the name() rather than a hard coded string 'diLepton' as the array index
	// std::string role = patCompCand.name();
	plots["diLepton"]["p"  ]->Fill(patCompCand.p());
	plots["diLepton"]["pt" ]->Fill(patCompCand.pt());
	plots["diLepton"]["eta"]->Fill(patCompCand.eta());
	plots["diLepton"]["phi"]->Fill(patCompCand.phi());
	plots["diLepton"]["charge"]->Fill(patCompCand.charge());
	plots["diLepton"]["mass"]->Fill(patCompCand.mass());
      }

      // now fill the daughter histograms
      // here we can use the CompositeCandidate role() method.  WHy not the name() method above?
      std::vector<std::string> role_collection = patCompCand.roles();
      //LogDebug("ZP2M")<<role_collection.size();
      //foreach(const std::string & str, role_collection) LogDebug("ZP2M")<<str;
      //LogDebug("ZP2M")<<"nDaughter "<<patCompCand.numberOfDaughters();
      for (size_t i = 0; i < patCompCand.numberOfDaughters(); ++i) {
	const reco::CandidateBaseRef dau = patCompCand.daughter(i)->masterClone();
	const pat::Muon* mu = toConcretePtr<pat::Muon>(dau);
	std::string role = role_collection[i];

	// do we want the selector to choose leptons or dileptons?
        if (!selector_(*mu)) continue;

	// basic kinematics
        plots[role]["p"  ]->Fill(mu->p());
        plots[role]["pt" ]->Fill(mu->pt());
        plots[role]["eta"]->Fill(mu->eta());
        plots[role]["phi"]->Fill(mu->phi());
        plots[role]["charge"]->Fill(mu->charge());

        if (mu->innerTrack().isNonnull()) {
            plots[role]["pixelHits"  ]->Fill(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
            plots[role]["pixelLayers"]->Fill(mu->innerTrack()->hitPattern().pixelLayersWithMeasurement());
            plots[role]["trackerHits"]->Fill(mu->innerTrack()->hitPattern().numberOfValidHits());
            plots[role]["trackerLostHitsMiddle"]->Fill(mu->innerTrack()->hitPattern().numberOfLostHits());
            plots[role]["trackerLostHitsInner"]->Fill(mu->innerTrack()->trackerExpectedHitsInner().numberOfLostHits());
            plots[role]["trackerLostHitsOuter"]->Fill(mu->innerTrack()->trackerExpectedHitsOuter().numberOfLostHits());
            plots[role]["trackerChi2n"]->Fill(mu->innerTrack()->normalizedChi2());

            if (!vertices->empty() && !vertices->front().isFake()) {
                const reco::Vertex &vtx = vertices->front();
                plots[role]["dxyCoarse"]->Fill(mu->innerTrack()->dxy(vtx.position()));
                plots[role]["dzCoarse"]->Fill(mu->innerTrack()->dz(vtx.position()));
                plots[role]["dxyFine"]->Fill(mu->innerTrack()->dxy(vtx.position()));
                plots[role]["dzFine"]->Fill(mu->innerTrack()->dz(vtx.position()));
            }
        }
        if (mu->outerTrack().isNonnull()) {
	    plots[role]["pSta"  ]->Fill(mu->outerTrack()->p());
	    plots[role]["ptSta" ]->Fill(mu->outerTrack()->pt());
	    plots[role]["etaSta"]->Fill(mu->outerTrack()->eta());
	    plots[role]["phiSta"]->Fill(mu->outerTrack()->phi());
	  
            plots[role]["muonHits"]->Fill(mu->outerTrack()->numberOfValidHits());
            plots[role]["muonChi2n"]->Fill(mu->outerTrack()->normalizedChi2());
            if ( ( mu->outerTrack()->extra().isAvailable()   ) && 
                 ( mu->outerTrack()->recHitsSize() > 0       ) &&
                 ( mu->outerTrack()->recHit(0).isAvailable() )     ) {
	        plots[role]["muonBadHits"]->Fill(mu->outerTrack()->recHitsSize() - mu->outerTrack()->numberOfValidHits());
	        plots[role]["muonStationsValid"]->Fill(mu->outerTrack()->hitPattern().muonStationsWithValidHits());
                plots[role]["muonStationsAny"  ]->Fill(mu->outerTrack()->hitPattern().muonStationsWithAnyHits());
                float abseta = std::abs(mu->outerTrack()->eta());
                if (abseta <= 1.2) {
		    plots[role]["muonStationsDTValid"]->Fill(mu->outerTrack()->hitPattern().dtStationsWithValidHits());
		    plots[role]["muonStationsDTAny"  ]->Fill(mu->outerTrack()->hitPattern().dtStationsWithAnyHits());
                } 
                if (abseta <= 1.6) {
		    plots[role]["muonStationsRPCValid"]->Fill(mu->outerTrack()->hitPattern().rpcStationsWithValidHits());
		    plots[role]["muonStationsRPCAny"  ]->Fill(mu->outerTrack()->hitPattern().rpcStationsWithAnyHits());
                } 
                if (abseta >= 0.8) {
		    plots[role]["muonStationsCSCValid"]->Fill(mu->outerTrack()->hitPattern().cscStationsWithValidHits());
		    plots[role]["muonStationsCSCAny"  ]->Fill(mu->outerTrack()->hitPattern().cscStationsWithAnyHits());
                }
            }
        }
        if (mu->globalTrack().isNonnull()) {
            plots[role]["globalHits"]->Fill(mu->globalTrack()->numberOfValidHits());
            plots[role]["globalMuonHits"]->Fill(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
            plots[role]["globalChi2n"]->Fill(mu->globalTrack()->normalizedChi2());
        }

        if (mu->isIsolationValid()) {
            plots[role]["trackIso05"]->Fill(mu->isolationR05().sumPt);
            plots[role][ "ecalIso05"]->Fill(mu->isolationR05().emEt);
            plots[role][ "hcalIso05"]->Fill(mu->isolationR05().hadEt);
            plots[role]["trackIso03"]->Fill(mu->isolationR03().sumPt);
            plots[role][ "ecalIso03"]->Fill(mu->isolationR03().emEt);
            plots[role][ "hcalIso03"]->Fill(mu->isolationR03().hadEt);
            plots[role][ "combRelIso03"]->Fill( (mu->isolationR03().sumPt + mu->isolationR03().emEt + mu->isolationR03().hadEt) / mu->pt() );
            plots[role][ "combRelIso05"]->Fill( (mu->isolationR05().sumPt + mu->isolationR05().emEt + mu->isolationR05().hadEt) / mu->pt() );
        }
        
        if (mu->isMatchesValid()) {
            plots[role]["segmentMatchesArb"    ]->Fill(mu->numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
            plots[role]["segmentMatchesNoArb"  ]->Fill(mu->numberOfMatches(reco::Muon::SegmentArbitration));
            plots[role]["segmentMatchesFailArb"]->Fill(mu->numberOfMatches(reco::Muon::SegmentArbitration) - mu->numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
            plots[role]["segmentCompatArb"     ]->Fill(muon::segmentCompatibility(*mu, reco::Muon::SegmentAndTrackArbitration));
            plots[role]["segmentCompatNoArb"   ]->Fill(muon::segmentCompatibility(*mu, reco::Muon::SegmentArbitration));
        }

        if (mu->isCaloCompatibilityValid()) {
            plots[role]["caloCompat"]->Fill(mu->caloCompatibility());
        }

       }
	
    }
	
}

void Zprime2muAnalysisPlots::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{
    if (luminosity != 0) {
        edm::Handle<edm::MergeableCounter> mc;
        iLumi.getByLabel(normalization_, mc);
        luminosity->Fill(0.5, double(mc->value));
        // set the correct uncertainty from counting statistics
        luminosity->SetBinError(1, sqrt(luminosity->GetBinContent(1)));
    }
}

DEFINE_FWK_MODULE(Zprime2muAnalysisPlots);







