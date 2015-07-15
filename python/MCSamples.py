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
    ##sample('zpsi5000',  'Z\'_{#psi} (5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM',      97890,  48, 1.,  0.0369,  k_factor=1.3, is_zprime=True),
    #sample('zpsi2250',  'Z\'_{#psi} (2.25 TeV) #rightarrow #mu^{+}#mu^{-}', '/RelValZpMM_13/CMSSW_7_4_0-MCRUN2_74_V7_gensim_740pre7-v1/GEN-SIM-RECO',      9000,  48, 1.,  0.0369,  k_factor=1.3, is_zprime=True),
    sample('dy50',   'DY50',   '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM',   19917018, 4 , 1., 6025.2,    k_factor=1.),
           ##sample('dy50',   'DY50',   '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM',   19917018, 37 , 1., 6104,    k_factor=1.027),5740 nic, 6025.2 nnlo
    #sample('dy120',  'DY120',  '/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM',   60756, 432, 1., 37,      k_factor=1.049),6025.2
    #  sample('dy1400', 'DY1400', '/DYJetsToEEMuMu_M-1400To2300_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 23204, 433, 1., 0.00255, k_factor=1.049),
    #  sample('dy200',  'DY200',  '/DYJetsToEEMuMu_M-200To400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',   26430, 434, 1., 5.26,    k_factor=1.049),
    #    sample('dy2300', 'DY2300', '/DYJetsToEEMuMu_M-2300To3500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM', 37358, 435, 1., 0.000165, k_factor=1.049),
    #    sample('dy3500', 'DY3500', '/DYJetsToEEMuMu_M-3500To4500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 23297, 436, 1., 0.00000758, k_factor=1.049),
    #    sample('dy400',  'DY400',  '/DYJetsToEEMuMu_M-400To800_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',   48205, 437, 1., 0.458,   k_factor=1.049),
    #    sample('dy4500', 'DY4500', '/DYJetsToEEMuMu_M-4500To6000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 25925, 438, 1., 0.00000796, k_factor=1.049),
    #    sample('dy6000', 'DY6000', '/DYJetsToEEMuMu_M-6000To7500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 28710, 439, 1., 0.000000022, k_factor=1.049),
    #    sample('dy7500', 'DY7500', '/DYJetsToEEMuMu_M-7500To8500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 30151, 431, 1., 0.000000000384, k_factor=1.049),
    #    sample('dy800',  'DY800',  '/DYJetsToEEMuMu_M-800To1400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',  24003, 3  , 1., 0.0314, k_factor=1.049),
    #    sample('dy8500', 'DY8500', '/DYJetsToEEMuMu_M-8500To9500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM', 32583, 36 , 1., 0.0000000000182, k_factor=1.049),
    #    sample('dy9500', 'DY9500', '/DYJetsToEEMuMu_M-9500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM',       37002, 8  , 1., 0.000000000000573, k_factor=1.049),
    sample('ttbar','t#bar{t}', '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM',                  4901956, 2 , 1., 832, k_factor=1.),
    #    sample('tW'   ,'tW',       '/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',    986100,  1, 1., 35.6, k_factor=1.),
    #    sample('tbarW','tbarW',    '/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 971800, 12, 1., 35.6, k_factor=1.),
    sample('ww',    'WW',      '/WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM',           499926,   5, 1., 12.59 , k_factor=1. ),
           #sample('ww',    'WW',      '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM',           697656,   4, 1., 63.21),
    sample('wz',    'WZ',      '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM',           996920, 6, 1., 66.1, k_factor=1.), #inclusive
    #    sample('zz'   ,'ZZ',       '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM',                  1958600, 6, 1., 15.4, k_factor=1.),
    sample('zz'   ,'ZZ',       '/ZZTo4L_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM',                  6621404, 9, 1., 0.157, k_factor=1.),
    sample('wjets','W+jets',   '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM',            22463202, 3, 1., 61500, k_factor=1.),
           ##sample('inclmu15', 'QCD',  '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM', 4767935, 801, 1., 867000000, k_factor=1.),
    ]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    #if '_c' in sample.name:
    #if 'Zprime' in sample.dataset:
#    sample.ana_dataset = '/%s/federica-%s-f646da20575c2cb2b2eda7b4413fb91e/USER'  % (sample.dataset.split('/')[1], sample.name)
    sample.ana_dataset = '/%s/rradogna-datamc_%s-c4b4ec8fa143ea00cec443e9d0afb38f/USER'  % (sample.dataset.split('/')[1], sample.name)

        #else:
            #sample.ana_dataset = '/%s/federica-%s-02dba98b5abbcd2765544ae02b3dcc74/USER'  % (sample.dataset.split('/')[1], sample.name) # this is actually wrong


ttbar.ana_dataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_ttbar-1e36332d8badf10b79a5027340f46eb1/USER'
wjets.ana_dataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_wjets-1e36332d8badf10b79a5027340f46eb1/USER'


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
