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
# Single-top cross sections are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma8TeV
# K factor for Drell-Yan samples is the ratio of the NNLO to POWHEG cross sections for M > 20 GeV bin, 1915/1871=1.024
samples = [
    sample('zmumu',     '#gamma/Z #rightarrow #mu^{+}#mu^{-}',              '/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',      3293740, 432, 0.05, 1915.),
#    sample('dy120',     'DY120',                                            '/DYToMuMu_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',       99984, 433, 0.05, 11.89,    k_factor=1.024),
#    sample('dy200',     'DY200',                                            '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',       99990, 434, 0.05, 1.485,    k_factor=1.024),
#    sample('dy500',     'DY500',                                            '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',       99992, 435, 0.05, 0.04415,  k_factor=1.024),
#    sample('dy800',     'DY800',                                            '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',       99984,   3, 0.05, 0.005491, k_factor=1.024),
#    sample('dy1000',    'DY1000',                                           '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',      99989,  36, 0.05, 0.001796, k_factor=1.024),
#    sample('dy1500',    'DY1500',                                           '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',      99992,   8, 0.05, 1.705E-4, k_factor=1.024),
#    sample('dy2000',    'DY2000',                                           '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',      99974,  37, 0.05, 2.208E-5, k_factor=1.024),
    sample('dy120_c1',  'DY120_C1',                                         '/DYToMuMu_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       99984, 433, 0.05, 11.89,    k_factor=1.024),
    sample('dy200_c1',  'DY200_C1',                                         '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       99990, 434, 0.05, 1.485,    k_factor=1.024),
    sample('dy500_c1',  'DY500_C1',                                         '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       99992, 435, 0.05, 0.04415,  k_factor=1.024),
    sample('dy800_c1',  'DY800_C1',                                         '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       99984,   3, 0.05, 0.005491, k_factor=1.024),
    sample('dy1000_c1', 'DY1000_C1',                                        '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',      99989,  36, 0.05, 0.001796, k_factor=1.024),
    sample('dy1500_c1', 'DY1500_C1',                                        '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',      99992,   8, 0.05, 1.705E-4, k_factor=1.024),
    sample('dy2000_c1', 'DY2000_C1',                                        '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',      99974,  37, 0.05, 2.208E-5, k_factor=1.024),
#    sample('dy120_c2',  'DY120_C2',                                         '/DYToMuMu_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',       99984, 433, 0.05, 11.89,    k_factor=1.024),
#    sample('dy200_c2',  'DY200_C2',                                         '/DYToMuMu_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',       99990, 434, 0.05, 1.485,    k_factor=1.024),
#    sample('dy500_c2',  'DY500_C2',                                         '/DYToMuMu_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',       99992, 435, 0.05, 0.04415,  k_factor=1.024),
#    sample('dy800_c2',  'DY800_C2',                                         '/DYToMuMu_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',       99984,   3, 0.05, 0.005491, k_factor=1.024),
#    sample('dy1000_c2', 'DY1000_C2',                                        '/DYToMuMu_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',      99989,  36, 0.05, 0.001796, k_factor=1.024),
#    sample('dy1500_c2', 'DY1500_C2',                                        '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',      99992,   8, 0.05, 1.705E-4, k_factor=1.024),
#    sample('dy2000_c2', 'DY2000_C2',                                        '/DYToMuMu_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM',      99974,  37, 0.05, 2.208E-5, k_factor=1.024),
#    sample('ttbar',     't#bar{t}',                                         '/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM', 6923750,   2, 0.067, 234.),
    sample('ttbar_powheg','t#bar{t}',                                        '/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM', 21675970,   2, 0.067, 234.),
    sample('tW',        'tW',                                               '/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',           497658,   1, 0.069, 11.1),
    sample('tbarW',     'tbarW',                                            '/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',        493460,  12, 0.069, 11.1),
    sample('ww',        'WW',                                               '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',                     10000431,   4, 0.035, 54.8),
    sample('wz',        'WZ',                                               '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',                     10000283,  30, 0.038, 33.2),
    sample('zz',        'ZZ',                                               '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',                      9799908,   6, 0.025,  8.1),
    sample('ztautau',   'Z #rightarrow #tau^{+}#tau^{-}',                   '/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',    3295238,  46, 0.05, 1915.),
    sample('wjets',     'W+jets',                                           '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',           18393090,   3, 0.05, 36257.),
    sample('inclmu15',  'QCD',                                              '/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',     21484602, 801, 0.1,  3.64E8 * 3.7E-4),
#    sample('zpsi750_c1',   'Z\'_{#psi} (0.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',        25040,  48, 0.05,  0.14,    k_factor=1.3, is_zprime=True),
#    sample('zpsi1000_c1',  'Z\'_{#psi} (1 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25040,  48, 0.05,  0.0369,  k_factor=1.3, is_zprime=True),
#    sample('zpsi1250_c1',  'Z\'_{#psi} (1.25 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-1250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25344,  48, 0.05,  0.0129,  k_factor=1.3, is_zprime=True),
#    sample('zpsi1500_c1',  'Z\'_{#psi} (1.5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimePSIToMuMu_M-1500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25344,  48, 0.05,  0.00433, k_factor=1.3, is_zprime=True),
#    sample('zpsi1750_c1',  'Z\'_{#psi} (1.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-1750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25272,  48, 0.05,  0.00172, k_factor=1.3, is_zprime=True),
#    sample('zpsi2000_c1',  'Z\'_{#psi} (2 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-2000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25092,  48, 0.05,  6.88E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2250_c1',  'Z\'_{#psi} (2.25 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-2250_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25104,  48, 0.05,  2.93E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2500_c1',  'Z\'_{#psi} (2.5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimePSIToMuMu_M-2500_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25344,  48, 0.05,  1.27E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2750_c1',  'Z\'_{#psi} (2.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-2750_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25376,  48, 0.05,  5.55E-5, k_factor=1.3, is_zprime=True),
#    sample('zpsi3000_c1',  'Z\'_{#psi} (3 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-3000_TuneZ2star_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C1-v1/AODSIM',       25040,  48, 0.05,  2.5E-5,  k_factor=1.3, is_zprime=True),
]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    if '_c' in sample.name:
        sample.ana_dataset = '/%s/slava-datamc_%s-7cd4d04801ad7f47970af9f536392613/USER' % (sample.dataset.split('/')[1], sample.name)
    else:
        sample.ana_dataset = '/%s/slava-datamc_%s-ecac376f8fa7ccc229aaa06d757d785a/USER' % (sample.dataset.split('/')[1], sample.name)

ttbar_powheg.ana_dataset  = '/TT_CT10_TuneZ2star_8TeV-powheg-tauola/slava-datamc_ttbar_powheg-7cd4d04801ad7f47970af9f536392613/USER'

#big_warn('nothing')
#big_warn('subtracting 232314 from ww.nevents because 3 tupling jobs got stuck')
#ww.nevents -= 232314

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
