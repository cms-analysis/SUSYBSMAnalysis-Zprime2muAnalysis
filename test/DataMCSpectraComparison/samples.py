#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import files_from_dbs

class sample:
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='REDIGI38X', ana_dataset=None):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
        self.syst_frac = syst_frac
        self.cross_section = cross_section
        self.k_factor = k_factor
        self.filenames_ = filenames
        self.scheduler = scheduler
        self.hlt_process_name = hlt_process_name
        self.ana_dataset_ = ana_dataset

    @property
    def partial_weight(self):
        return self.cross_section / float(self.nevents) * self.k_factor # the total weight is partial_weight * integrated_luminosity

    @property
    def ana_dataset(self):
        if self.ana_dataset_ is not None:
            return self.ana_dataset_
        publish_log_fn = 'crab/publish_logs/publish.crab_datamc_%s' % self.name
        ad = [x.strip().replace('=== dataset ', '') for x in open(publish_log_fn).readlines() if x.startswith('=== dataset')]
        assert(len(ad) == 1)
        self.ana_dataset_ = ad[0]
        return ad[0]

    @property
    def filenames(self):
        # Return a list of filenames for running the histogrammer not
        # using crab.
        if self.filenames_ is not None:
            return self.filenames_
        return files_from_dbs(self.ana_dataset, ana02=True)

    def __getitem__(self, key):
        return getattr(self, key)

    def _dump(self, redump_existing=False):
        dst = os.path.join('/uscmst1b_scratch/lpc1/3DayLifetime/tucker', self.name)
        os.system('mkdir ' + dst)
        for fn in self.filenames:
            print fn
            if redump_existing or not os.path.isfile(os.path.join(dst, os.path.basename(fn))):
                os.system('dccp ~%s %s/' % (fn,dst))
            
# https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries for
# xsecs (all below in pb)
samples = [
    sample('zmumu',        'Z #rightarrow #mu^{+}#mu^{-}',             '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-START38_V12-v1/GEN-SIM-RECO',          1998931, 432, 0.096, 1631 - 0.97*1.3, hlt_process_name='HLT'),
    sample('dy200',        'DY200',                                    '/DYToMuMu_M-200_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                             55000, 433, 0.1,   0.97  - 0.027, k_factor=1.3, hlt_process_name='HLT'),
    sample('dy500',        'DY500',                                    '/DYToMuMu_M-500_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                             55000, 434, 0.1,   0.027 - 0.003, k_factor=1.3, hlt_process_name='HLT'),
    sample('dy800',        'DY800',                                    '/DYToMuMu_M-800_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                             55000, 435, 0.1,   0.003,          k_factor=1.3, hlt_process_name='HLT'),
    sample('ttbar',        't#bar{t}',                                 '/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO',                    1167759,   2, 0.15,  157,  hlt_process_name='HLT'),
    sample('singletop_tW', 'tW',                                       '/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/Fall10-START38_V12-v2/GEN-SIM-RECO',                494961,   1, 0.075, 10.6, hlt_process_name='HLT'),
    sample('ww',           'WW',                                       '/WWtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2061760,   4, 0.035, 43,  scheduler='glite'),
    sample('wz',           'WZ',                                       '/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2194752,   5, 0.038, 18,  scheduler='glite'),
    sample('zz',           'ZZ',                                       '/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2113368,   6, 0.025, 5.9, scheduler='glite'),
    sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',           '/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO', 1994719,  46, 0.096, 1631, hlt_process_name='HLT'),
    sample('wjets',        'W+jets',                                   '/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               15168266,   3, 0.086, 1.04e4, hlt_process_name='HLT'),
    sample('inclmu15',     'QCD (MuRich, muon p_{T} > 15 GeV)',        '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',        29504866, 801, 0.1,   0.0002855 * 296600000),
    sample('zssm750',    'Z\'_{SSM}(750) #rightarrow #mu^{+}#mu^{-}',  '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                      55000,  40, 0.1,   0.355, k_factor=1.3, hlt_process_name='HLT'),
]
samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn
big_warn('modifying ttbar and singletop_tW nevents down by 25k since you lost one job from each')
ttbar.nevents -= 25000
singletop_tW.nevents -= 25000
big_warn('also modifying wz down by 25k and and wjets down by 50k since 1/2 files of the tuple are stuck')
wz.nevents -= 25000
wjets.nevents -= 25000

print 'samples are:'
print ' '.join(s.name for s in samples)
