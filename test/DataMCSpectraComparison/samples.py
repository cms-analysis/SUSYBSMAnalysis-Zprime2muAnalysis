#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import files_from_dbs

class sample:
    def __init__(self, name, nice_name, dataset, nevents, color, cross_section, k_factor=1, filenames=None, scheduler='condor', is_35x=True, hlt_process_name='REDIGI', ana_dataset=None):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
        self.cross_section = cross_section
        self.k_factor = k_factor
        self.filenames_ = filenames
        self.scheduler = scheduler
        self.is_35x = is_35x
        self.hlt_process_name = hlt_process_name
        self.partial_weight = cross_section / float(nevents) * k_factor # the total weight is partial_weight * integrated_luminosity
        self.ana_dataset_ = ana_dataset

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


# https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries for
# xsecs (all below in pb)
samples = [
    sample('zmumu',        'Z #rightarrow #mu^{+}#mu^{-}',         '/Zmumu_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v2/GEN-SIM-RECO',          1768457,   7, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('ttbar',        't#bar{t}',                             '/TTbarJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',              1483404,   2, 152),
    sample('singletop_tW', 'Single t (tW)',                        '/SingleTop_tWChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',     466437,   1, 10.6),
    sample('ww',           'WW',                                   '/WW/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               122980,   4, 43),
    sample('wz',           'WZ',                                   '/WZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               118120,   5, 18),
    sample('zz',           'ZZ',                                   '/ZZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               145368,   6, 5.9),
    sample('wjets',        'W+jets',                               '/WJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                 10068895,   3, 2.8e4), # is this xsec right?
    sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',       '/Ztautau_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v1/GEN-SIM-RECO',        1269404,  46, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('inclmu15',     'QCD',                                          '/InclusiveMu15/Summer10-START37_V5_S09-v1/GEN-SIM-RECO',                    5369781, 801, 2.969e8*0.00037, is_35x=False, hlt_process_name='REDIGI37X'),
    sample('zssm750',      "Z_{SSM} (750) #rightarrow #mu^{+}#mu^{-}",     '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO',    18932,  40, 0.355, hlt_process_name='HLT'),
    sample('qcd15',        'QCD (15 < #hat{p_{T}} < 20 GeV)',    '/QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',    2884915, 801, 0.00254 * 579200000, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd20',        'QCD (20 < #hat{p_{T}} < 30 GeV)',    '/QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',   11461085, 801, 0.00518 * 236300000, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd30',        'QCD (30 < #hat{p_{T}} < 50 GeV)',    '/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',   11431864, 801, 0.0109  *  53070000, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd50',        'QCD (50 < #hat{p_{T}} < 80 GeV)',    '/QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',   10748755, 801, 0.02274 *   6351000, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd80',        'QCD (80 < #hat{p_{T}} < 120 GeV)',   '/QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',   3191979, 801, 0.037   *    785100, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd120',       'QCD (120 < #hat{p_{T}} < 150 GeV)',  '/QCD_Pt-120to150_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',   998503, 801, 0.04777 *     92950, hlt_process_name='REDIGI38X', is_35x=False),
    sample('qcd150',       'QCD (#hat{p_{T}} > 150 GeV)',        '/QCD_Pt-150_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',       1022541, 801, 0.05964 *     47580, hlt_process_name='REDIGI38X', is_35x=False),
]
samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn
big_warn('modifying wjets nevents (by just a small amount, though)')
wjets.nevents = 10067072

print 'samples are:'
print ' '.join(s.name for s in samples)
