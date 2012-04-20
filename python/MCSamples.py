#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn, files_from_dbs
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import dataset_from_publish_log

class sample(object):
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='HLT', ana_dataset=None, is_zprime=False):
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
        self.ana_dataset = ana_dataset
        self.is_zprime = is_zprime

    @property
    def partial_weight(self):
        return self.cross_section / float(self.nevents) * self.k_factor # the total weight is partial_weight * integrated_luminosity

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

class tupleonlysample(sample):
    def __init__(self, name, dataset, scheduler='condor', hlt_process_name='HLT'):
        super(tupleonlysample, self).__init__(name, 'dummy', dataset, 1, 1, 1, 1, scheduler=scheduler, hlt_process_name=hlt_process_name)

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV for xsecs (all below in pb)
# commented out samples are not available for Summer12/8 TeV
samples = [
    sample('zmumu',        '#gamma/Z #rightarrow #mu^{+}#mu^{-}',                '/DYToMuMu_M_20_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START50_V15-v1/AODSIM',               1963296, 432, 0.05, 1915 - 0), # The subtraction here and below is because we drop the events with generated M in ranges in the other datasets, e.g. drop here M > 120.
#   sample('dy120',        'DY120',                                              '/DYToMuMu_M-120_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',             54550, 433, 0.05, 7.9    - 0.97,   k_factor=1.3),
#   sample('dy200',        'DY200',                                              '/DYToMuMu_M-200_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',             55000, 433, 0.05, 0.97   - 0.027,  k_factor=1.3),
#   sample('dy500',        'DY500',                                              '/DYToMuMu_M-500_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',             55000, 434, 0.05, 0.027  - 0.0031, k_factor=1.3),
#   sample('dy800',        'DY800',                                              '/DYToMuMu_M-800_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',             55000, 435, 0.05, 0.0031 - 9.7e-4, k_factor=1.3),
#   sample('dy1000',       'DY1000',                                             '/DYToMuMu_M-1000_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',            55000, 435, 0.05, 9.7e-4,          k_factor=1.3),
    sample('ttbar',        't#bar{t}',                                           '/TTJets_TuneZ2star_8TeV-madgraph-tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',              1203373,   2, 0.15, 225),
#   sample('tW',           'tW',                                                 '/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM',            814390,   1, 0.03,   7.9),
#   sample('tbarW',        'tbarW',                                              '/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM',         809984,   1, 0.03,   7.9),
    sample('ww',           'WW',                                                 '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',                  10000431,   4, 0.035, 57, scheduler='glite'),
    sample('wz',           'WZ',                                                 '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',                  10000283,  30, 0.038, 32, scheduler='glite'),
    sample('zz',           'ZZ',                                                 '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',                   9799908,   6, 0.025,  8, scheduler='glite'),
#   sample('ztautau',      'Z #rightarrow #tau^{+}#tau^{-}',                     '/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola/Summer11-PU_S3_START42_V11-v2/AODSIM',          2032536,  46, 0.05, 1631),
#   sample('wjets',        'W+jets',                                             '/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM',             81352581,   3, 0.05, 3.1e4),
#   sample('inclmu15',     'QCD',                                                '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM',      25080241, 801, 0.1,   2.855e-4 * 2.966e8),
#   sample('zssm1000',     'Z\'_{SSM} (1 TeV) #rightarrow #mu^{+}#mu^{-}',       '/ZprimeSSMToMuMu_M-1000_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM',            20328,  38, 0.05,  0.089, k_factor=1.3, is_zprime=True),
#   sample('zpsi1000',     'Z\'_{#psi} (1 TeV) #rightarrow #mu^{+}#mu^{-}',      '/ZprimePSIToMuMu_M-1000_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM',            21020,  48, 0.05,  0.089, k_factor=1.3, is_zprime=True),
]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    sample.ana_dataset = '/%s/tucker-datamc_%s-caca636ea661546409f4073c061b3e20/USER' % (sample.dataset.split('/')[1], sample.name)

#big_warn('nothing')

__all__ = ['samples'] + [s.name for s in samples]


if __name__ == '__main__':
    if False:
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            parents = dbsparents(s.dataset)
            for parent in parents:
                for line in os.popen('dbss rel %s' % parent):
                    if 'CMSSW' in line:
                        print parent, line,
            print

    if False:
        import os
        from dbstools import dbsparents
        for s in [ww,wz,zz]:
            print s.dataset
            parents = dbsparents(s.dataset)
            print parents
            os.system('dbsconfig %s > %s' % (parents[-1], s.name))

        os.system('dbss nevents %s' % x.replace('RECO','RAW'))
        os.system('dbss nevents %s' % x)

    if False:
        import os
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            def fuf(y):
                x = os.popen(y).read()
                for line in x.split('\n'):
                    try:
                        print int(line)
                    except ValueError:
                        pass
            fuf('dbss nevents %s' % s.dataset)
            fuf('dbss nevents %s' % s.dataset.replace('AODSIM','GEN-SIM-RECO'))

    if False:
        for s in samples:
            print s.name
            os.system('grep "total events" ~/nobackup/crab_dirs/384p3/publish_logs/publish.crab_datamc_%s' % s.name)
            os.system('grep "total events" ~/nobackup/crab_dirs/413p2/publish_logs/publish.crab_datamc_%s' % s.name)
            print

    if False:
        os.system('mkdir ~/scratch/wjets')
        for fn in wjets.filenames:
            assert fn.startswith('/store')
            fn = '/pnfs/cms/WAX/11' + fn
            cmd = 'dccp %s ~/scratch/wjets/' % fn
            print cmd
            os.system(cmd)

    if False:
        for s in samples:
            print s.name
            os.system('dbss site %s' % s.dataset)
            print

    if False:
        for s in samples:
            if s.ana_dataset is None:
                continue
            c = []
            for line in os.popen('dbss ana02 find file.numevents where dataset=%s' % s.ana_dataset):
                try:
                    n = int(line)
                except ValueError:
                    continue
                c.append(n)
            c.sort()
            print s.name, c
