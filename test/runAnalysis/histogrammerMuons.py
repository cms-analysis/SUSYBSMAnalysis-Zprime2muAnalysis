#!/usr/bin/env python
Electrons = False

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
####################################
####################################
####################################

HistosFromPAT.usekFactor = True #### Set TRUE to use K Factor #####
HistosFromPAT.useTTBarWeight = True #### Set TRUE to use NNPDF Weights for ttbar #####

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
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff as OurSelection2018






# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
	'Our2017'  : OurSelection2017,
	'Our2018'  : OurSelection2018,
	}
if year == 2016:
	cuts = {
		'Our2016'  : OurSelection2016,
		'Our2017'  : OurSelection2017,
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
    # present.
    path_list.append(process.egmGsfElectronIDSequence)
	    
   
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut

    leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
    if year == 2016 and (isMC or "03Feb" in sampleName or "23Sep" in sampleName or "Prompt" in sampleName):
	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')

    if len(trigger_filters)>0 and (cut_name=='Our2017' or cut_name=='Simple' or cut_name == 'Our2018'):
    	leptons.trigger_filters = trigger_filters
	leptons.trigger_path_names = trigger_path_names
        leptons.trigger_path_full_names = trigger_path_full_names
	leptons.prescaled_trigger_filters = prescaled_trigger_filters
	leptons.prescaled_trigger_path_names = prescaled_trigger_path_names
    if len(trigger_filters)>0 and year == 2016:
    	leptons.trigger_filters = trigger_filters2016
	leptons.trigger_path_names = trigger_path_names2016
        leptons.trigger_path_full_names = trigger_path_full_names2016
	leptons.prescaled_trigger_filters = prescaled_trigger_filters
	leptons.prescaled_trigger_path_names = prescaled_trigger_path_names


#    if isMC:
#	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')
    if  Electrons:
	    if cut_name == 'EmuVeto':
		    leptons.electron_muon_veto_dR = 0.1

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
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dimuons.clone(src = cms.InputTag(allname))
	if len(trigger_filters) >  0 and (cut_name=='Our2017' or cut_name=='Our2016' or cut_name=='Simple'):
		if year == 2016:
			alldil.tight_cut = trigger_match_2016
		else:
			alldil.tight_cut = trigger_match_2018
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
      
	histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name), year = cms.int32(year))
        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
	if not isMC:
		delattr(getattr(process,name + 'Histos'),'hardInteraction')	

        path_list.append(alldil * dil * histos)

	if 'ConLR' in sampleName or 'DesLR' in sampleName or 'ConRL' in sampleName or 'DesRL' in sampleName:
		L = 10000	
		if '16TeV' in sampleName:
			L = 16000
		if '100kTeV' in sampleName:
			L = 100000000
		if '1TeV' in sampleName:
			L = 1000
		if '22TeV' in sampleName:
			L = 22000
		if '24TeV' in sampleName:
			L = 24000
		if '28TeV' in sampleName:
			L = 28000
		if '32TeV' in sampleName:
			L = 32000
		if '34TeV' in sampleName:
			L = 34000
		if '40TeV' in sampleName:
			L = 40000
		histos.lrWeightProducer.Lambda = L	
		histos.lrWeightProducer.calculate = True
		histos.lrWeightProducer.doingElectrons = False
		if "RL" in sampleName:
			histos.lrWeightProducer.doingLR = False
		if 'Con' in sampleName:
			histos.lrWeightProducer.interference = -1
		else:	
			histos.lrWeightProducer.interference = 1
	
    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
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


if addNTuples:

	process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
					   dimu_src = cms.InputTag('Our2018MuonsPlusMuonsMinus'),
						met_src = cms.InputTag("slimmedMETs"),
						jet_src = cms.InputTag("slimmedJets"),
					   beamspot_src = cms.InputTag('offlineBeamSpot'),
					   vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
	# 								TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
								TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
					   genEventInfo = cms.untracked.InputTag('generator'),
					   metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter")),
					   doElectrons = cms.bool(False),
					   )
	if isMC:
		process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
		obj = process.prunedMCLeptons
		obj.src = cms.InputTag('prunedGenParticles')

		from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction_MiniAOD
		process.SimpleNtupler.hardInteraction = hardInteraction_MiniAOD
		if year == 2016:
			process.SimpleNtupler.dimu_src = cms.InputTag('Our2016MuonsPlusMuonsMinus')
			if hasattr(process, 'pathOur2016'):
				process.pathOur2016 *=obj * process.SimpleNtupler 
		else:
			if hasattr(process, 'pathOur2018'):
				process.pathOur2018 *=obj * process.SimpleNtupler 
	else:
		if year == 2016:
			process.SimpleNtupler.dimu_src = cms.InputTag('Our2016MuonsPlusMuonsMinus')
			if hasattr(process, 'pathOur2016'):
				process.pathOur2016 *= process.SimpleNtupler 
		else:
			if hasattr(process, 'pathOur2018'):
				process.pathOur2018 *=  process.SimpleNtupler 

if isMC:
	switch_reco_process_name(process, "PAT") # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)


print process.dumpPython()
