#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =['/store/relval/CMSSW_10_6_0_pre4/SingleMuon/MINIAOD/106X_dataRun2_tracker2017_v1_RelVal_sigMu2017F-v1/20000/9D181EE9-0CD9-4E40-BF15-42C5860C4C68.root']

process.maxEvents.input = -1
isMC = True
addNTuples = False
year = 2017
sampleName = 'dummy'
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

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
Electrons = False

from SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionAtZ_cfi import ResolutionAtZ
ResolutionAtZ.leptonsFromDileptons = True
ResolutionAtZ.doQoverP = True

from SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi import ResolutionUsingMC_MiniAOD as ResolutionUsingMC
ResolutionUsingMC.leptonsFromDileptons = True
ResolutionUsingMC.doQoverP = True


####################################
####################################
####################################


####################################
####################################
####################################


####################################
####################################
####################################

ZSkim = False #### Set TRUE to skim dy50to120 with a Z pt < 100 GeV #####

####################################
####################################
####################################

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2017_cff as OurSelection2017






# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	#('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	#('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	#('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	#('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
	'Our2016'  : OurSelection2016,
	'Our2017'  : OurSelection2017,
	#'OurNoIso' : OurSelectionDec2012,
	#'Simple'   : OurSelection2017, # The selection cuts in the module listed here are ignored below.
	#'OurMuPrescaledNew'  : OurSelectionNew,
	#'OurMuPrescaled2012' : OurSelectionDec2012
	}

tracks = ["Inner","Global","Outer","TPFMS","DYT","Picky","TunePNew"]
# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for trackType in tracks:
	for cut_name, Selection in cuts.iteritems():
		# Keep track of modules to put in the path for this set of cuts.
	    path_list = []

	    # Clone the LeptonProducer to make leptons with the set of cuts
	    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
	    # just marks each muon with a userInt "cutFor" that is 0 if it
	    # passes the cuts, and non-0 otherwise; it does not actually drop
	    # any of the muons. The cutFor flag actually gets ignored by the
	    # LooseTightPairSelector in use for all the cuts above, at
	    # present.
	    path_list.append(process.egmGsfElectronIDSequence)
		    
            leptons_name = cut_name + 'Leptons' + trackType
	    if cut_name == 'Simple':
		muon_cuts = ''
#	    elif 'MuPrescaled' in cut_name:
	    muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > 20')
#	    else:
#		muon_cuts = Selection.loose_cut

	    leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
	    leptons.muon_track_for_momentum = cms.string(trackType)
	#    if isMC:
	#	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')
	    if year == 2016 and isMC:
		leptons.trigger_summary = cms.InputTag('selectedPatTrigger')


	    if  Electrons:
		    if cut_name == 'EmuVeto':
			    leptons.electron_muon_veto_dR = 0.1
	    if len(trigger_filters)>0 and (cut_name=='Our2017' or cut_name=='Simple'):
		leptons.trigger_filters = trigger_filters
		leptons.trigger_path_names = trigger_path_names
		leptons.prescaled_trigger_filters = prescaled_trigger_filters
		leptons.prescaled_trigger_path_names = prescaled_trigger_path_names

	    # Keep using old TuneP for past selections
	#     if 'Dec2012' not in Selection.__file__:
	#         leptons.muon_track_for_momentum = cms.string('TunePNew')
	    setattr(process, leptons_name, leptons)
	    path_list.append(leptons)

	    # Make all the combinations of dileptons we defined above.
	    for dil_name, dil_decay, dil_cut in dils:
		# For the EmuVeto path, we only care about e-mu events.
		if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
		    continue

		# For the MuPrescaled paths, we don't care about e-mu events.
		if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
		    continue

		# Unique names for the modules: allname for the allDileptons,
		# and name for dileptons.
		name = cut_name + dil_name + trackType	
		allname = 'all' + name

		alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
		if 'AllSigns' in dil_name:
		    alldil.checkCharge = cms.bool(False)

		dil = Selection.dimuons.clone(src = cms.InputTag(allname))
		#if len(trigger_filters) >  0 and (cut_name=='Our2017' or cut_name=='Simple'):
		alldil.tight_cut = ""
		alldil.loose_cut = muon_cuts
		# Implement the differences to the selections; currently, as
		# in Zprime2muCombiner, the cuts in loose_cut and
		# tight_cut are the ones actually used to drop leptons, and
		# not the ones passed into the LeptonProducer to set cutFor above.
		if cut_name == 'Simple':
		    alldil.electron_cut_mask = cms.uint32(0)
		    #alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
		    alldil.loose_cut = 'isGlobalMuon && pt > 20.'
		    alldil.tight_cut = ''
		    dil.max_candidates = 100
		    dil.sort_by_pt = True
		    dil.do_remove_overlap = False
		    if hasattr(dil, 'back_to_back_cos_angle_min'):
			delattr(dil, 'back_to_back_cos_angle_min')
		    if hasattr(dil, 'vertex_chi2_max'):
			delattr(dil, 'vertex_chi2_max')
		    if hasattr(dil, 'dpt_over_pt_max'):
			delattr(dil, 'dpt_over_pt_max')
		elif 'MuPrescaled' in cut_name:
			alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
			assert alldil.tight_cut == trigger_match
			if len(prescaled_trigger_filters)>0:
				alldil.tight_cut = prescaled_trigger_match_2018
			else:
				alldil.tight_cut = prescaled_trigger_match
	    # Histos now just needs to know which leptons and dileptons to use.
	      
		res    = ResolutionAtZ.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))
		resMC  = ResolutionUsingMC.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))


		# Add all these modules to the process and the path list.
		setattr(process, allname, alldil)
		setattr(process, name, dil)
		setattr(process, name + 'Resolution', res)
		if isMC:
			setattr(process, name + 'ResolutionMC', resMC)
			path_list.append(alldil * dil * res * resMC)
		else:	
			path_list.append(alldil * dil * res)


	    # Finally, make the path for this set of cuts.
	    pathname = 'path' + cut_name + trackType
	    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
	    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
	    pobj = process.EventCounter * process.dileptonPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)



	    if 'VBTF' not in cut_name and cut_name != 'Simple':
		process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
		for dataFilter in goodDataFiltersMiniAOD:
			#setattr(process,dataFilter 
			pobj = dataFilter * pobj


	    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
		pobj = process.PrescaleToCommonMiniAOD * pobj ####### Now it seams that there are no prescaled path ########
	    path = cms.Path(pobj)
	    setattr(process, pathname, path)




if isMC:
	switch_reco_process_name(process, "PAT") # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
print process.dumpPython()
