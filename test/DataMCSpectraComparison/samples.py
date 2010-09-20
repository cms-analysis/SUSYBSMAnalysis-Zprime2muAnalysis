#!/usr/bin/env python

import os

class sample:
    @classmethod
    def psethash(cls, dataset):
        # These need to be specified with whatever crab made for the
        # published dataset name.
        if 'Zmumu' in dataset or 'Ztautau' in dataset:
            return 'ea1b4401edd0c9e8af9e80917519ee4e'
        return 'f1606b2f2aa07e70082d78e786896133'

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
        self.ana_dataset = '/%s/%s-datamc_%s-%s/USER' % (dataset.split('/')[1], os.environ['USER'], name, self.psethash(dataset)) if ana_dataset is None else ana_dataset

    @property
    def filenames(self):
        # Return a list of filenames for running the histogrammer not
        # using crab.
        if self.filenames_ is not None:
            return self.filenames_
        # could use DBSAPI but this is easier
        cmd = 'dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query "find file where dataset=%s"' % self.ana_dataset
        return [y.strip('\n') for y in os.popen(cmd).readlines() if '.root' in y]

# https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries for
# xsecs (all below in pb)
samples = [
    sample('zmumu',        'Z #rightarrow #mu^{+}#mu^{-}',         '/Zmumu_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v2/GEN-SIM-RECO',          1768457,   7, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('ttbar',        't#bar{t}',                             '/TTbarJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',              1483404,   2, 152),
    sample('singletop_tW', 'Single t (tW)',                        '/SingleTop_tWChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',     466437,   1, 10.6),
#   sample('singletop_s',  'Single t (s-channel)',                 '/SingleTop_sChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',      412055,  30, 4.6*0.32442),
#   sample('singletop_t',  'Single t (t-channel)',                 '/SingleTop_tChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',      528593,  40, 63*0.32442),
    sample('wjets',        'W+jets',                               '/WJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                 10068895,   3, 2.8e4), # is this xsec right?
    sample('ww',           'WW',                                   '/WW/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               122980,   4, 43),
    sample('wz',           'WZ',                                   '/WZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               118120,   5, 18),
    sample('zz',           'ZZ',                                   '/ZZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',                               145368,   6, 5.9),
    sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',       '/Ztautau_M20_CTEQ66-powheg/Summer10-START36_V9_S09-v1/GEN-SIM-RECO',        1269404,  46, 1686, is_35x=False, hlt_process_name='REDIGI36X'),
    sample('qcd100',       'QCD (100 < #hat{p}_{T} < 250 GeV/c)',  '/QCD_Pt100to250-madgraph/Spring10-START3X_V26_S09-v2/GEN-SIM-RECO',        10875132, 801, 7e6),    # LO xsec!
    sample('qcd250',       'QCD (250 < #hat{p}_{T} < 500 GeV/c)',  '/QCD_Pt250to500-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',         4913036, 802, 1.71e5), # LO xsec!
    sample('qcd500',       'QCD (500 < #hat{p}_{T} < 1000 GeV/c)', '/QCD_Pt500to1000-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',        4234762, 803, 5e3),    # LO xsec!
    sample('qcd1000',      'QCD (#hat{p}_{T} > 1000 GeV/c)',       '/QCD_Pt1000toInf-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO',        1661261, 804, 80),     # LO xsec!
]
samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name

def warn(s):
    x = '#' * len(s)
    print x
    print x
    print x
    print s
    print x
    print x
    print x

warn('modifying qcd100,250 and wjets nevents (by just a small amount, though))
wjets.nevents = 10067072
qcd100.nevents = 10862502
qcd250.nevents -= 25000

print 'samples are:'
for i,s in enumerate(samples):
    print '%2i:%20s' % (i, s.name)
