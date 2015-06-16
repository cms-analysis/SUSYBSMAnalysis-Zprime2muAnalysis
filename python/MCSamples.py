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
# check cross sections (dibosons)
samples = [
           sample('zpsi5000',  'Z\'_{#psi} (5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM',      97890,  2, 1.,  5.48E-5,  k_factor=1., is_zprime=True),
           sample('dy50',   'DY50',   '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',   2829164, 37 , 1., 5740,    k_factor=1.),
           sample('dy120',  'DY120',  '/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM',   60756, 432, 1., 37,      k_factor=1.),
           sample('dy1400', 'DY1400', '/DYJetsToEEMuMu_M-1400To2300_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 23204, 433, 1., 7.58E-6 , k_factor=1.),
           sample('dy200',  'DY200',  '/DYJetsToEEMuMu_M-200To400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',   26430, 434, 1., 2.55E-3,    k_factor=1.),
           sample('dy2300', 'DY2300', '/DYJetsToEEMuMu_M-2300To3500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM', 37358, 435, 1., 4.58E-1, k_factor=1.),
           sample('dy3500', 'DY3500', '/DYJetsToEEMuMu_M-3500To4500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 23297, 436, 1., 7.96E-6, k_factor=1.),
           sample('dy400',  'DY400',  '/DYJetsToEEMuMu_M-400To800_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',   48205, 432, 1., 5.26,   k_factor=1.),
           sample('dy4500', 'DY4500', '/DYJetsToEEMuMu_M-4500To6000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 25925, 433, 1., 2.20E-8, k_factor=1.),
           sample('dy6000', 'DY6000', '/DYJetsToEEMuMu_M-6000To7500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 28710, 434, 1., 3.84E-10, k_factor=1.),
           sample('dy7500', 'DY7500', '/DYJetsToEEMuMu_M-7500To8500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 30151, 435, 1., 3.14E-2, k_factor=1.),
           sample('dy800',  'DY800',  '/DYJetsToEEMuMu_M-800To1400_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',  24003, 436, 1., 1.65E-4, k_factor=1.),
           sample('dy8500', 'DY8500', '/DYJetsToEEMuMu_M-8500To9500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM', 32583, 36 , 1., 1.82E-11, k_factor=1.),
           sample('dy9500', 'DY9500', '/DYJetsToEEMuMu_M-9500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v2/AODSIM',       37002, 8  , 1., 5.73E-13, k_factor=1.),
           sample('ttbar','t#bar{t}', '/TT_Tune4C_13TeV-pythia8-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM',                  2991609, 2 , 1., 832, k_factor=1.),
           sample('tW'   ,'tW',       '/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',    986100,  1, 1., 35.6, k_factor=1.),
           sample('tbarW','tbarW',    '/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM', 971800, 12, 1., 35.6, k_factor=1.),
           #    sample('ww',    'WW',      '/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',              899900,   4, 1., 1.08E+1, k_factor=1.),
           sample('wz'   ,'WZ',       '/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',           237484, 30, 1., 1.634, k_factor=1.),
           sample('zz'   ,'ZZ',       '/ZZTo4L_Tune4C_13TeV-powheg-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',                  1958600, 5, 1., 15.4, k_factor=1.),
           sample('wjets','W+jets',   '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM',            10017930, 2, 1., 61500, k_factor=1.),
           sample('inclmu15', 'QCD',  '/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v3/AODSIM', 4767935, 801, 1., 8.67E+8, k_factor=1.),
           ]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    #if '_c' in sample.name:
        #if 'ZprimePSI' in sample.dataset:
            #sample.ana_dataset = '/%s/valuev-datamc_%s-7cd4d04801ad7f47970af9f536392613/USER' % (sample.dataset.split('/')[1], sample.name)
        #else:
            #sample.ana_dataset = '/%s/slava-datamc_%s-7cd4d04801ad7f47970af9f536392613/USER' % (sample.dataset.split('/')[1], sample.name)
    #else:
    
    if 'Zprime' in sample.dataset:
        #sample.ana_dataset = '/%s/federica-%s-bcef9cff6c54eb003cec7431c0bde17e/USER'  % (sample.dataset.split('/')[1], sample.name)
        sample.ana_dataset = '/%s/rradogna-datamc_%s-a8881ceec144e0dfafbb7486d1b7f8e6/USER' % (sample.dataset.split('/')[1], sample.name)
    else:
        sample.ana_dataset = '/%s/federica-%s-02dba98b5abbcd2765544ae02b3dcc74/USER'  % (sample.dataset.split('/')[1], sample.name) # this is actually wrong


#big_warn('nothing')
#big_warn('subtracting 232314 from ww.nevents because 3 tupling jobs got stuck')
#ww.nevents -= 232314

#zz.ana_dataset = '/ZZTo4L_Tune4C_13TeV-powheg-pythia8/federica-ZZ_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#wz.ana_dataset = '/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/federica-WZ_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#wjets.ana_dataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/federica-WJets_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#tbarW.ana_dataset = '/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/federica-Tbar_tW_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#tW.ana_dataset = '/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/federica-T_tW_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#ttbar.ana_dataset = '/TT_Tune4C_13TeV-pythia8-tauola/federica-TT_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#inclmu15.ana_dataset = '/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/federica-QCD_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy50.ana_dataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/federica-DY_M-50_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy9500.ana_dataset = '/DYJetsToEEMuMu_M-9500_13TeV-madgraph/federica-DY_M-9500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy8500.ana_dataset = '/DYJetsToEEMuMu_M-8500To9500_13TeV-madgraph/federica-DY_M-8500To9500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy800.ana_dataset = '/DYJetsToEEMuMu_M-800To1400_13TeV-madgraph/federica-DY_M-800To1400_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy7500.ana_dataset = '/DYJetsToEEMuMu_M-7500To8500_13TeV-madgraph/federica-DY_M-7500To8500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy6000.ana_dataset = '/DYJetsToEEMuMu_M-6000To7500_13TeV-madgraph/federica-DY_M-6000To7500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy4500.ana_dataset = '/DYJetsToEEMuMu_M-4500To6000_13TeV-madgraph/federica-DY_M-4500To6000_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy400.ana_dataset = '/DYJetsToEEMuMu_M-400To800_13TeV-madgraph/federica-DY_M-400To800_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy3500.ana_dataset = '/DYJetsToEEMuMu_M-3500To4500_13TeV-madgraph/federica-DY_M-3500To4500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy2300.ana_dataset = '/DYJetsToEEMuMu_M-2300To3500_13TeV-madgraph/federica-DY_M-2300To3500_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy200.ana_dataset = '/DYJetsToEEMuMu_M-200To400_13TeV-madgraph/federica-DY_M-200To400_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy1400.ana_dataset = '/DYJetsToEEMuMu_M-1400To2300_13TeV-madgraph/federica-DY_M-1400To2300_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'
#dy120.ana_dataset = '/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/federica-DY_M-120To200_Phys14_PU20BX25-02dba98b5abbcd2765544ae02b3dcc74/USER'

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
