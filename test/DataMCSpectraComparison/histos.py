#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT

import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection_cff as OurSelection

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
common_dil_cut = '' #daughter(0).momentum.Dot(daughter(1).momentum())/daughter(0).momentum().Mag()/daughter(1).momentum().Mag() > 0.02'
dils = [
    ('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    ('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
    ('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
    ('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
    ('ElectronsPlusElectronsMinus',  '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@-', 'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    ('ElectronsPlusElectronsPlus',   '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@+', 'daughter(0).pdgId() + daughter(1).pdgId() == -22'),
    ('ElectronsMinusElectronsMinus', '%(leptons_name)s:electrons@- %(leptons_name)s:electrons@-', 'daughter(0).pdgId() + daughter(1).pdgId() == 22'),
    ('ElectronsSameSign',            '%(leptons_name)s:electrons@- %(leptons_name)s:electrons@-', ''),
    ('MuonsPlusElectronsMinus',      '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == -2'),
    ('MuonsMinusElectronsPlus',      '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == 2'),
    ('MuonsPlusElectronsPlus',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == -24'),
    ('MuonsMinusElectronsMinus',     '%(leptons_name)s:muons@- %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == 24'),
    ('MuonsElectronsOppSign',        '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     ''),
    ('MuonsElectronsSameSign',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
    ]

# Define groups of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add hltFilter
# somewhere below.
cuts = ['VBTF', 'Our', 'OurNoIso', 'OurIso3', 'OurRelIso015', 'OurRelIso006', 'OurNoPx']

for cut_name in cuts:
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.
    leptons_name = cut_name + 'Leptons'
    muon_cuts = VBTFSelection.loose_cut if 'VBTF' in cut_name else OurSelection.loose_cut
    leptons = process.leptons.clone(muon_cuts = muon_cuts)
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        name = cut_name + dil_name
        allname = 'all' + name

        if common_dil_cut and dil_cut:
            dil_cut = common_dil_cut + ' && (%s)' % dil_cut

        if 'VBTF' in cut_name:
            alldil, dil = VBTFSelection.allDimuons, VBTFSelection.dimuons
        else:
            alldil, dil = OurSelection.allDimuons, OurSelection.dimuons
            
        alldil = alldil.clone(decay = dil_decay % locals(), cut = dil_cut)
        dil = dil.clone(src = cms.InputTag(allname))

        if 'NoIso' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt < 10', '')
        elif 'Iso3' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt < 10', ' && isolationR03.sumPt < 3')
        elif 'RelIso015' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt < 10', ' && isolationR03.sumPt / innerTrack.pt < 0.15')
        elif 'RelIso006' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt < 10', ' && isolationR03.sumPt / innerTrack.pt < 0.06')
        elif 'NoPx' in cut_name:
            alldil.tight_cut = alldil.tight_cut.value().replace(' && innerTrack.hitPattern.numberOfValidPixelHits >= 1', '')
                    
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)

    # Finally, make the path for this set of cuts. Don't use hltFilter
    # here, but rely on the selection to take care of it -- easy way
    # to handle the changing trigger names.
    pathname = 'path' + cut_name
    path = cms.Path(process.goodDataFilter * process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list))
    setattr(process, pathname, path)

if 'olddata' in sys.argv:
    process.source.fileNames = ['file:work/daata/jul15.root', 'file:work/daata/prompt.root']
    process.TFileService.fileName = 'ana_datamc_data.root'
    process.GlobalTag.globaltag = 'GR10_P_V7::All'
elif 'data' in sys.argv:
    process.source.fileNames = ['file:crab/crab_datamc_Run2010A/res/merged.root', 'file:crab/crab_datamc_promptB_all/res/merged.root']
    process.TFileService.fileName = 'ana_datamc_data.root'    
    process.GlobalTag.globaltag = 'GR10_P_V10::All'
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*eval(open('Run2010AB.cmssw').read()))

def ntuplify(process, hlt_process_name='HLT'):
    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler', hlt_src = cms.InputTag('TriggerResults', '', hlt_process_name), dimu_src = cms.InputTag('OurMuonsPlusMuonsMinus'))
    process.SimpleNtuplerSS = process.SimpleNtupler.clone(dimu_src = cms.InputTag('OurMuonsSameSign'))
    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('OurMuonsElectronsOppSign'))
    process.pathOur *= process.SimpleNtupler * process.SimpleNtuplerSS * process.SimpleNtuplerEmu

    process.SimpleNtuplerVBTF = process.SimpleNtupler.clone(dimu_src = cms.InputTag('VBTFMuonsPlusMuonsMinus'))
    process.pathVBTF *= process.SimpleNtuplerVBTF

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')
    process.PrintEvent = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('OurMuonsPlusMuonsMinus'))
    process.PrintEventSS = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurMuonsSameSign'))
    process.PrintEventEmu = process.PrintEvent.clone(dilepton_src = cms.InputTag('OurMuonsElectronsOppSign'))
    process.pathOur *= process.PrintEvent * process.PrintEventSS * process.PrintEventEmu

    process.PrintEventVBTF = process.PrintEvent.clone(dilepton_src = cms.InputTag('VBTFMuonsPlusMuonsMinus'))
    process.pathVBTF *= process.PrintEventVBTF
    
if 'data' in sys.argv or 'olddata' in sys.argv:
    printify(process)
    ntuplify(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(ana_dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = histos_crab.py
get_edm_output = 1
total_number_of_events = -1
events_per_job = 20000

[USER]
ui_working_dir = crab/crab_ana_datamc_%(name)s
return_data = 1
'''

    just_testing = 'testing' in sys.argv

    from samples import samples
    for sample in samples:
        print sample.name
        
        new_py = open('histos.py').read()
        new_py += "\nprocess.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', '%(hlt_process_name)s')\n" % sample
        new_py += "\nntuplify(process, hlt_process_name='%(hlt_process_name)s')\n" % sample
            
        open('histos_crab.py', 'wt').write(new_py)

        open('crab.cfg', 'wt').write(crab_cfg % sample)
        if not just_testing:
            os.system('crab -create -submit all')
        
    if not just_testing:
        os.system('rm crab.cfg histos_crab.py histos_crab.pyc')
