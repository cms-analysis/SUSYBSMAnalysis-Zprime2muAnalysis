#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import files_from_dbs

class sample:
    def __init__(self, name, nice_name, dataset, nevents, color, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='REDIGI38X', ana_dataset=None):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
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


# https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries for
# xsecs (all below in pb)
samples = [
    sample('zmumu',        'Z #rightarrow #mu^{+}#mu^{-}',             '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-START38_V12-v1/GEN-SIM-RECO',          1998931,   7, 1631, hlt_process_name='HLT'),
    sample('ttbar',        't#bar{t}',                                 '/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO',                    1167759,   2, 152,  hlt_process_name='HLT'),
    sample('singletop_tW', 'tW',                                       '/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/Fall10-START38_V12-v2/GEN-SIM-RECO',                494961,   1, 10.6, hlt_process_name='HLT'),
    sample('ww',           'WW',                                       '/WWtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2061760,   4, 43, scheduler='glite'),
    sample('wz',           'WZ',                                       '/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2194752,   5, 18, scheduler='glite'),
    sample('zz',           'ZZ',                                       '/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO',               2113368,   6, 5.9, scheduler='glite'),
    sample('wmunu',        'W #rightarrow #mu#nu',                     '/WToMuNu_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                           5330940,   3, 7899), # do not include the filter eff listed in the file because it wasn't run on this production -- to be checked
    sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',           '/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO', 1994719,  46, 1631, hlt_process_name='HLT'),
    sample('inclmu15',     'QCD (MuRich, muon p_{T} > 15 GeV)',        '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',        29504866, 801, 0.0002855 * 296600000),
#   sample('qcd15',        'QCD (15 < #hat{p_{T}} < 20 GeV)',          '/QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',       2884915, 801, 0.00254 * 579200000),
#   sample('qcd20',        'QCD (20 < #hat{p_{T}} < 30 GeV)',          '/QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',      11461085, 801, 0.00518 * 236300000),
#   sample('qcd30',        'QCD (30 < #hat{p_{T}} < 50 GeV)',          '/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',      11431864, 801, 0.0109  *  53070000),
#   sample('qcd50',        'QCD (50 < #hat{p_{T}} < 80 GeV)',          '/QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',      10748755, 801, 0.02274 *   6351000),
#   sample('qcd80',        'QCD (80 < #hat{p_{T}} < 120 GeV)',         '/QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',      3191979, 801, 0.037   *    785100),
#   sample('qcd120',       'QCD (120 < #hat{p_{T}} < 150 GeV)',        '/QCD_Pt-120to150_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',      998503, 801, 0.04777 *     92950),
#   sample('qcd150',       'QCD (#hat{p_{T}} > 150 GeV)',              '/QCD_Pt-150_MuPt5Enriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',          1022541, 801, 0.05964 *     47580),
#   sample('qcdem20',      'QCD (EMRich, 20 < #hat{p_{T}} < 30 GeV)',  '/QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',         37169939, 802, 0.0104  * 236000000),
#   sample('qcdem30',      'QCD (EMRich, 30 < #hat{p_{T}} < 80 GeV)',  '/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',         71845473, 803, 0.065   *  59480000),
#   sample('qcdem80',      'QCD (EMRich, 80 < #hat{p_{T}} < 170 GeV)', '/QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',         8073559, 804, 0.155   *    900000),
#   sample('qcdbce20',     'QCD (BCtoE, 20 < #hat{p_{T}} < 30 GeV)',   '/QCD_Pt-20to30_BCtoE_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',               2243439, 805, 0.00056 * 236000000),
#   sample('qcdbce30',     'QCD (BCtoE, 30 < #hat{p_{T}} < 80 GeV)',   '/QCD_Pt-30to80_BCtoE_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',               1995502, 806, 0.0023  *  59480000),
#   sample('qcdbce80',     'QCD (BCtoE, 80 < #hat{p_{T}} < 170 GeV)',  '/QCD_Pt-80to170_BCtoE_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',              1043390, 807, 0.0104  *    900000),
    sample('zssm750',      'Z_{SSM} (750) #rightarrow #mu^{+}#mu^{-}', '/ZprimeSSMToMuMu_M-750_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO',                      55000,  40, 0.355, hlt_process_name='HLT'),
]
samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn
big_warn('modifying ttbar and singletop_tW nevents down by 25k since you lost one job from each')
ttbar.nevents -= 25000
singletop_tW.nevents -= 25000

print 'samples are:'
print ' '.join(s.name for s in samples)








