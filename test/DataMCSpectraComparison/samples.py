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

# https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries for
# xsecs (all below in pb)
samples = [
    sample('zmumu',        'Z #rightarrow #mu^{+}#mu^{-}',         '/Zmumu_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v2/GEN-SIM-RECO',          1768457,   7, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('ttbar',        't#bar{t}',                             '/TTbarJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',              1483404,   2, 152),
    sample('singletop_tW', 'Single t (tW)',                        '/SingleTop_tWChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',     466437,   1, 10.6),
#   sample('singletop_s',  'Single t (s-channel)',                 '/SingleTop_sChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',      412055,  30, 4.6*0.32442),
#   sample('singletop_t',  'Single t (t-channel)',                 '/SingleTop_tChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',      528593,  40, 63*0.32442),
    sample('ww',           'WW',                                   '/WW/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               122980,   4, 43),
    sample('wz',           'WZ',                                   '/WZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               118120,   5, 18),
    sample('zz',           'ZZ',                                   '/ZZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               145368,   6, 5.9),
    sample('wjets',        'W+jets',                               '/WJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                 10068895,   3, 2.8e4), # is this xsec right?
    sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',       '/Ztautau_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v1/GEN-SIM-RECO',        1269404,  46, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('qcd100',       'QCD (100 < #hat{p}_{T} < 250 GeV)',    '/QCD_Pt100to250-madgraph/Spring10-START3X_V26_S09-v2/GEN-SIM-RECO',        10875132, 801, 7e6),    # LO xsec!
    sample('qcd250',       'QCD (250 < #hat{p}_{T} < 500 GeV)',    '/QCD_Pt250to500-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',         4913036, 802, 1.71e5), # LO xsec!
    sample('qcd500',       'QCD (500 < #hat{p}_{T} < 1000 GeV)',   '/QCD_Pt500to1000-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',        4234762, 803, 5e3),    # LO xsec!
    sample('qcd1000',      'QCD (#hat{p}_{T} > 1000 GeV)',         '/QCD_Pt1000toInf-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',        1661261, 804, 80),     # LO xsec!
]
samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn
big_warn('modifying qcd100,250 and wjets nevents (by just a small amount, though)')
wjets.nevents = 10067072
qcd100.nevents = 10862502
qcd250.nevents -= 25000

print 'samples are:'
print ' '.join(s.name for s in samples)
