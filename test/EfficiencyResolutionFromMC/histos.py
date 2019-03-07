#!/usr/bin/env python

################################################################################

miniAOD = True
check_prescaled_path = False
# User beware: note do_all_track_fits=True can potentially cause memory problems on CRAB
# Remove unnecessary EDAnalyzers or cut sets from the path
do_all_track_fits = False 
ex = ''

################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples

#process.maxEvents.input = 1000
process.source.fileNames = [
            '/store/mc/RunIIAutumn18MiniAOD/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/2BBCC514-3E29-1D40-BCA5-B0F0D09FD08C.root',
                            ]
process.options.wantSummary = False

process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v12'

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True

from SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi import ResolutionUsingMC
ResolutionUsingMC.leptonsFromDileptons = cms.bool(True)
ResolutionUsingMC.doQoverP = cms.bool(True)
ResolutionUsingMC.hardInteraction.src = cms.InputTag('prunedGenParticles')

from SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi import EfficiencyFromMCMini as EfficiencyFromMC
EfficiencyFromMC.use_resonance_mass = cms.bool(True)
EfficiencyFromMC.use_resonance_mass_denom = cms.bool(True)

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process

if do_all_track_fits:
    tracks = [('muons',''),('global','Global'), ('inner','Inner'), ('tpfms','TPFMS'), ('picky','Picky'), ('dyt','DYT')]
else:
    tracks = [('muons','')]

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold, trigger_filters, trigger_path_names, prescaled_trigger_filters, prescaled_trigger_path_names, prescaled_trigger_match_2018, trigger_match_2018

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff as OurSelection2018
cuts = {
    #'Our2012'  : OurSelectionDec2012,
    #'Our2016'  : OurSelection2016,
    'Our2018'  : OurSelection2018,
    'Our2018AtZ' : OurSelection2018,
    # Vertex-constrained mass not guaratneed with Simple selection
    #'Simple'   : OurSelection2018, # The selection cuts in the module listed here are ignored below.
    }
dils = [
        ('OppSign','%(leptons_name)s:%(track)s@+ %(leptons_name)s:%(track)s@-','daughter(0).pdgId() + daughter(1).pdgId() == 0')
        ]


process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)
process.HardInteractionFilterRes.hardInteraction.src = cms.InputTag('prunedGenParticles')

for cut_name, Selection in cuts.iteritems():

    path_list = []

    leptons_name = cut_name+'Leptons'
    if miniAOD: path_list.append(process.egmGsfElectronIDSequence)
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'AtZ' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut

    leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
    leptons.trigger_filters = trigger_filters
    leptons.trigger_path_names = trigger_path_names
    leptons.prescaled_trigger_filters = prescaled_trigger_filters
    leptons.prescaled_trigger_path_names = prescaled_trigger_path_names
    # put different tev muon reconstructors into event
    tracks_for_momentum = [track for track,name in tracks if track!="muons"] # muons is default TuneP
    leptons.muon_tracks_for_momentum = cms.vstring(tracks_for_momentum)

    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Loop on different tev muon track reconstructors
    for track,track_name in tracks:

        # Make all the combinations of dileptons we defined above.
        for dil_name, dil_decay, dil_cut in dils:

            # Unique names for the modules: allname for the allDileptons,
            # and name for dileptons.
            name = cut_name + track_name + dil_name
            allname = 'all' + name

            alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
            dil = Selection.dimuons.clone(src = cms.InputTag(allname))

            # Implement the differences to the selections; currently, as
            # in Zprime2muCombiner, the cuts in loose_cut and
            # tight_cut are the ones actually used to drop leptons, and
            # not the ones passed into the LeptonProducer to set cutFor above.
            if cut_name == 'Simple':
                alldil.electron_cut_mask = cms.uint32(0)
                alldil.loose_cut = 'isGlobalMuon && pt > 20.'
                alldil.tight_cut = ''
                dil.max_candidates = 100
                dil.sort_by_pt = True
                dil.do_remove_overlap = False
                dil.prefer_Z=False
                if hasattr(dil, 'back_to_back_cos_angle_min'):
                    delattr(dil, 'back_to_back_cos_angle_min')
                if hasattr(dil, 'vertex_chi2_max'):
                    delattr(dil, 'vertex_chi2_max')
                if hasattr(dil, 'dpt_over_pt_max'):
                    delattr(dil, 'dpt_over_pt_max')
            elif 'AtZ' in cut_name:
                alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
                alldil.tight_cut = prescaled_trigger_match_2018

            # Histos now just needs to know which leptons and dileptons to use.

            Histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name,'muons'),dilepton_src=cms.InputTag(name))

            ResolutionUsingMC.lepton_src = cms.InputTag(leptons_name,'muons')
            ResolutionUsingMC.dilepton_src = cms.InputTag(name)
            Resolution = ResolutionUsingMC.clone()
            ResolutionVertex = ResolutionUsingMC.clone(use_vertex_mass=cms.bool(True))

            if 'AtZ' in cut_name:
                EfficiencyFromMC.trigger_filters = prescaled_trigger_filters
                EfficiencyFromMC.trigger_path_names = prescaled_trigger_path_names
                EfficiencyFromMC.dimuon_src = cms.InputTag(name)
                Efficiency = EfficiencyFromMC.clone()
            else:
                EfficiencyFromMC.trigger_filters = trigger_filters
                EfficiencyFromMC.trigger_path_names = trigger_path_names
                EfficiencyFromMC.dimuon_src = cms.InputTag(name)
                Efficiency = EfficiencyFromMC.clone(trigger_filters=trigger_filters, trigger_path_names=trigger_path_names)

            setattr(process,allname,alldil)
            setattr(process,name,dil)
            setattr(process,name+'Histos',Histos)
            setattr(process,name+'Resolution',Resolution)
            setattr(process,name+'ResolutionVertex',ResolutionVertex)
            setattr(process,name+'Efficiency',Efficiency)
            path_list.append(alldil * dil * Histos * Resolution * ResolutionVertex * Efficiency)

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    pobj = process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)
    path = cms.Path(pobj)
    setattr(process, pathname, path)

#f = file('outfile', 'w')
#f.write(process.dumpPython())
#f.close()

###############################################################################
    
import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ana_effres_%(name)s%(extra)s'
config.General.workArea = 'crab'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'

config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
config.Data.outputDatasetTag = 'ana_effres_%(name)s'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dimuon/2018/effres'
config.Site.storageSite = 'T2_CH_CERN'
                          
'''
        
    # Only do DY MC samples and apply HardInteractionFilter
    FilterInfo = {
        #'name':        (  lo,    hi),
        'dy50to120':    (  50,   120),
        'dy120to200':   ( 120,   200),
        'dy200to400':   ( 200,   400),
        'dy400to800':   ( 400,   800),
        'dy800to1400':  ( 800,  1400),
        'dy1400to2300': (1400,  2300),
        'dy2300to3500': (2300,  3500),
        'dy3500to4500': (3500,  4500),
        'dy4500to6000': (4500,  6000),
        'dy6000toInf':  (6000, 10000),
        }
    dySamples = [sample for sample in samples if sample.name in FilterInfo.keys()]

    just_testing = 'testing' in sys.argv

    for sample in dySamples:

        name = sample.name
        dataset = sample.dataset
        lo = FilterInfo[name][0]
        hi = FilterInfo[name][1]
        print name,lo,hi,dataset
        extra = '_'+ex if ex!='' else ''

        open('crabConfig.py', 'wt').write(crab_cfg % locals())

        new_py = open('histos.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = %i\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = %i\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = %i\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = %i\n' % hi
        open('histos_crab.py', 'wt').write(new_py)
        
        if not just_testing:
            os.system('crab submit -c crabConfig.py')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
