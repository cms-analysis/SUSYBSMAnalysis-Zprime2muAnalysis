
#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =["/store/mc/RunIISummer16MiniAODv2/ADDGravToLL_LambdaT-10000_M-1700_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/66BD020A-31CA-E611-B77E-0CC47A6C183A.root"] 
#!/usr/bin/env python
#import FWCore.ParameterSet.Config as cms
#from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
#from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
#from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
#from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

#process.source.fileNames =[#'file:./pat.root'
#'/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root',
#'/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F4380BF8-CBCF-E611-8891-0CC47A546E5E.root'
#'/store/mc/RunIIFall17MiniAODv2/ZToEE_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/8A226AAF-AC43-E811-AEF0-0CC47A4D764C.root',
#'file:pickevents.root',
#'/store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/50000/000DCB8B-2ADD-E711-9100-008CFAF35AC0.root',
#'/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M1300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/DCA196D6-1B79-E711-BADA-90B11C282313.root'
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FEE03939-E778-E711-976D-008CFAF5550C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4E0B4D52-E678-E711-BFE6-0242AC1C0502.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/162AE54B-E678-E711-AFDF-48FD8EE73AC5.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CC2887A8-E778-E711-B955-3417EBE34C27.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F049644B-E678-E711-9036-0025905D1E08.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/20ED2453-E678-E711-9711-0CC47A4C8F30.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2E_M1300_CUETP8M1_Lam28TeVDesLR_13TeV_Pythia8_Corrected-v4/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/84357A56-E678-E711-8F3E-0CC47A4D764C.root",
#'/store/mc/PhaseIFall16MiniAOD/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/FlatPU28to62HcalNZSRAW_PhaseIFall16_exo52_90X_upgrade2017_realistic_v6_C1-v1/120000/304E419F-CC13-E711-93E9-FA163E0231A1.root',
#'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/12BD4CC4-0751-E811-BCA9-0090FAA58D84.root',
#'/store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam24TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/90F99FA7-E76B-E811-BC35-A0369FD0B362.root'
#]                          

process.maxEvents.input = -1
isMC = True
addNTuples = False
year = 2016
#process.GlobalTag.globaltag = '84X_mcRun2_asymptotic_2016_TranchelV_v6'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold, trigger_filters, trigger_path_names, prescaled_trigger_filters, prescaled_trigger_path_names, prescaled_trigger_match_2018, trigger_match_2018

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)
#!/usr/bin/env python
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor on DY. If used, the k factor will be applied to ALL samples submitted. #####

	
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
import SUSYBSMAnalysis.Zprime2muAnalysis.ElectronSelection_cff as ElectronSelection
import SUSYBSMAnalysis.Zprime2muAnalysis.ElectronSelection2016_cff as ElectronSelection2016






dils = [
	('ElectronsOppSign',        '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@-',     ''),
	('ElectronsSameSign',       '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@+',     ''),
	('ElectronsAllSigns',       '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@+',     ''),
	]

cuts = {
	'ElectronSelection'  : ElectronSelection,
	}
	
if year == 2016:
	cuts = {
	"ElectronSelection" : ElectronSelection2016,
	}
# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
	# Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present
    path_list.append(process.egmGsfElectronIDSequence)
	    
    leptons_name = cut_name + 'Leptons'
    leptons = process.leptonsMini.clone()
    if year == 2016 and isMC:
	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')

    setattr(process, leptons_name, leptons)
    path_list.append(leptons)
	
    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDielectrons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dielectrons.clone(src = cms.InputTag(allname))
	
	# Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        # Histos now just needs to know which leptons and dileptons to use.
	if isMC:
		if year == 2018:
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("mc_2017","data_2018"))
		elif year == 2017:	
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("mc_2017","data_2017"))
		else:	
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("mc_2016","data_2016"))
	else:	
		histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True))
	#if not isMC:
	#	delattr(histos,'hardInteraction')

	histos.hardInteraction.doingElectrons = True
        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
	if not isMC:
		#del histos.hardInteraction
		#histos.useMadgraphWeight  = False
		trig = Selection.dielectronHLT
		trigName = cut_name + "HLTFilter"
        	setattr(process, trigName, trig)
		delattr(getattr(process,name + 'Histos'),'hardInteraction')	
        	path_list.append(trig * alldil * dil * histos)
	else:
		alldil.loose_cut_ele = cms.string('et > 35 && abs(userFloat("etaSC")) < 2.5 && !(abs(userFloat("etaSC")) > 1.4442 && abs(userFloat("etaSC")) < 1.566)')	
		alldil.tight_cut_ele = cms.string("")	
		alldil.ele_match_l1 = cms.bool(False)	
        	path_list.append(alldil * dil * histos)
	
    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DielectronPreselector_cfi')
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
    pobj = process.EventCounter * process.dielectronPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)


    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
    for dataFilter in goodDataFiltersMiniAOD:
	#setattr(process,dataFilter 
	pobj = dataFilter * pobj

    path = cms.Path(pobj)
    setattr(process, pathname, path)


if addNTuples:
	  
	process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
					   dimu_src = cms.InputTag('ElectronSelectionElectronsAllSigns'),
					   met_src = cms.InputTag("slimmedMETs"),
					   jet_src = cms.InputTag("slimmedJets"),
					   beamspot_src = cms.InputTag('offlineBeamSpot'),
					   vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					   #TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
					   TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
					   genEventInfo = cms.untracked.InputTag('generator'),
					   metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter")),
					   doElectrons = cms.bool(True),
					   )


	if isMC:
		#process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
		#obj = process.prunedMCLeptons
		#obj.src = cms.InputTag('prunedGenParticles')
	 
		from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction_MiniAOD as hardInteraction
		hardInteraction.doingElectrons = True
		process.SimpleNtupler.hardInteraction = hardInteraction
		if hasattr(process, 'pathElectronSelection'):
			#process.pathElectronSelection *=obj * process.SimpleNtupler 
			process.pathElectronSelection *= process.SimpleNtupler 

	else:
		if hasattr(process, 'pathElectronSelection'):
			process.pathElectronSelection *= process.SimpleNtupler 


if isMC:
	switch_reco_process_name(process, "PAT") # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

process.DYGenMassFilter = cms.EDFilter('QScaleSelector',
                                       src = cms.InputTag('generator'),
                                       min_mass = cms.double(1700),
                                       max_mass = cms.double(7000)
                                       )
for path_name, path in process.paths.iteritems():
        getattr(process,path_name).insert(2,process.DYGenMassFilter)