#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn, files_from_dbs
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import dataset_from_publish_log

miniAOD = True

class sample(object):
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='HLT', ana_dataset=None, is_madgraph=False, is_zprime=False):
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
        self.is_madgraph = is_madgraph
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

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV for xsecs (all below in pb)
# Single-top cross sections are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
# K factor for Drell-Yan samples is the ratio of the NNLO to POWHEG cross sections for M > 20 GeV bin, 1915/1871=1.024
samples = [
    ####filippo's branch    
    #sample('dy50to120',   'DY50to120', '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 2967200, 209 , 1., 1975,   k_factor=1.006),#NLO xs and k-factor applied to reach NLO
    #sample('dy120to200',  'DY120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 99200, 210, 1., 19.32, k_factor=1.006),#mcm 19.32
    #sample('dy200to400',  'DY200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 211, 1., 2.731, k_factor=1.006),#mcm 2.731
    #sample('dy400to800',  'DY400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM', 100000, 212, 1., 0.241, k_factor=1.006),
    #sample('dy800to1400', 'DY800to1400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 72, 1., 0.01678, k_factor=1.006),
    #sample('dy1400to2300','DY1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 0.00139,    k_factor=1.006),
    #sample('dy2300to3500','DY2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 99200, 70 , 1., 0.00008948,    k_factor=1.006),
    #sample('dy3500to4500','DY3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 0.0000041,    k_factor=1.006),
    #sample('dy4500to6000','DY4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM', 100000, 70 , 1., 4.56E-7,    k_factor=1.006),

    sample('DY50to120Powheg',   'DY50to120',  '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM',  2976600, 209, 1., 1975, k_factor=1.006,hlt_process_name="HLT2"),#mcm 19.32
    sample('DY120to200Powheg',  'DY120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 100000, 210, 1., 19.32, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY200to400Powheg',  'DY200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 100000, 211, 1., 2.731, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY400to800Powheg',  'DY400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 99000, 212, 1., 0.241, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY800to1400Powheg', 'DY800to1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 96398, 72, 1., 0.01678, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY1400to2300Powheg','DY1400to2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 100000, 70, 1., 0.00139, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY2300to3500Powheg','DY2300to3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 100000, 70, 1., 0.00008948, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY3500to4500Powheg','DY3500to4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 99000, 70, 1., 0.0000041, k_factor=1.006,hlt_process_name="HLT2"),
    sample('DY4500to6000Powheg','DY4500to6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 99000, 70, 1., 4.56E-7, k_factor=1.006,hlt_process_name="HLT2"),  
    sample('DY6000toInfPowheg','DY6000toInf','/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 100000, 70, 1., 4.56E-7, k_factor=1.006,hlt_process_name="HLT2")    
    #sample('DY200to400Powheg',  'DY200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM', 100000, 211, 1., 2.731, k_factor=1.),#mcm 2.731
    #sample('DY400to800Powheg',  'DY400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM', 99000, 212, 1., 0.241, k_factor=1.,hlt_process_name="HLT2"),
    ]

samples.reverse()


#if miniAOD:

#else:
for sample in samples:
   exec '%s = sample' % sample.name
   if not miniAOD:
       sample.ana_dataset = '/%s/rradogna-datamc_%s-c4b4ec8fa143ea00cec443e9d0afb38f/USER'  % (sample.dataset.split('/')[1], sample.name)
   else:
       sample.ana_dataset = '/'+ sample.dataset.split('/')[1]+ '/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
# DY120to200Powheg.ana_dataset = 'ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
       print sample.ana_dataset
       #print  (sample.dataset.split('/AODSIM')[0]+ '/MINIAODSIM', sample.name)

            #sample.ana_dataset = '/%s/federica-%s-02dba98b5abbcd2765544ae02b3dcc74/USER'  % (sample.dataset.split('/')[1], sample.name) # this is actually wrong
#dy100to200.ana_dataset = '/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy100to200-1e36332d8badf10b79a5027340f46eb1/USER'
#dy200to400.ana_dataset = '/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy200to400-1e36332d8badf10b79a5027340f46eb1/USER'
#dy400to500.ana_dataset = '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy400to500-1e36332d8badf10b79a5027340f46eb1/USER'
#dy500to700.ana_dataset = '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy500to700-1e36332d8badf10b79a5027340f46eb1/USER'
#dy700to800.ana_dataset = '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy700to800-1e36332d8badf10b79a5027340f46eb1/USER'
#dy800to1000.ana_dataset = '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy800to1000-1e36332d8badf10b79a5027340f46eb1/USER'
#dy1000to1500.ana_dataset = '/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy1000to1500-1e36332d8badf10b79a5027340f46eb1/USER'
#dy1500to2000.ana_dataset = '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy1500to2000-1e36332d8badf10b79a5027340f46eb1/USER'
#dy2000to3000.ana_dataset = '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy2000to3000-1e36332d8badf10b79a5027340f46eb1/USER'

#zpsi5000.ana_dataset = '/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/rradogna-datamc_zpsi5000-1e36332d8badf10b79a5027340f46eb1/USER'
#dy50_startup.ana_dataset = '/RelValZMM_13/rradogna-datamc_dy50_startup-8b577364235a1c7c11f4fb31512a2917/USER'

#dy50to120.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/rradogna-datamc_dy50to120-1e36332d8badf10b79a5027340f46eb1/USER'
#DY120to200Powheg.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/alfloren-DY120to200Powheg-ea459820ba8ecaf0b251c44e2defe317/USER'
#DY200to400Powheg.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/alfloren-DY200to400Powheg-ea459820ba8ecaf0b251c44e2defe317/USER'
#DY400to800Powheg.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/alfloren-DY400to800Powheg-ea459820ba8ecaf0b251c44e2defe317/USER'
#DY800to1400Powheg.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/alfloren-DY800to1400Powheg-d361b004739dfc1dad40e50368455d7a/USER'
#dy1400to2300.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/rradogna-datamc_dy1400to2300-1e36332d8badf10b79a5027340f46eb1/USER'
#DY3500to4500Powheg.ana_dataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/alfloren-DY3500to4500Powheg-ea459820ba8ecaf0b251c44e2defe317/USER'

#ww_incl.ana_dataset = '/WW_TuneCUETP8M1_13TeV-pythia8/rradogna-datamc_ww_incl-1e36332d8badf10b79a5027340f46eb1/USER'
#zz_incl.ana_dataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/rradogna-datamc_zz_incl-1e36332d8badf10b79a5027340f46eb1/USER'
###wjets.ana_datasetOLD = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_wjets-1e36332d8badf10b79a5027340f46eb1/USER'
###ttbar.ana_dataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_ttbar-1e36332d8badf10b79a5027340f46eb1/USER'
#wjets.ana_dataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_wjets-d2059c7d2b57376da41472544da161fa/USER'
#ttbar_pow.ana_dataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/rradogna-datamc_ttbar_pow-1e36332d8badf10b79a5027340f46eb1/USER'
##ttbar_startup.ana_dataset = '/RelValTTbar_13/rradogna-datamc_ttbar_startup-8b577364235a1c7c11f4fb31512a2917/USER'
#tWantitop.ana_dataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/alfloren-tWantitop-728a04705e311faf7e2183c346d6b42c/USER'
#tWtop.ana_dataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/alfloren-tWtop-728a04705e311faf7e2183c346d6b42c/USER'


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
