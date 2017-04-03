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
    sample('dy50to120',   'DY50to120', '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 2977600, 209 , 1., 1975,   k_factor=1.),#NLO xs and k-factor applied to reach NLO
    sample('dy120to200',  'DY120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 210, 1., 19.32, k_factor=1.),#mcm 19.32
    sample('dy200to400',  'DY200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 211, 1., 2.731, k_factor=1.),#mcm 2.731
    sample('dy400to800',  'DY400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 98400, 212, 1., 0.241, k_factor=1.),
    sample('dy800to1400', 'DY800to1400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 72, 1., 0.01678, k_factor=1.),
    sample('dy1400to2300','DY1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 70 , 1., 0.00139,    k_factor=1.),
    sample('dy2300to3500','DY2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 70 , 1., 0.00008948,    k_factor=1.),
    sample('dy3500to4500','DY3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 70 , 1., 0.0000041,    k_factor=1.),
    sample('dy4500to6000','DY4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 100000, 70 , 1., 4.56E-7,    k_factor=1.),    

    
     sample('WZ', 'WZ', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 1000000, 98, 1., 47.13, k_factor=1.),#NLO from MCFM
     sample('WZ_ext', 'WZ_ext', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM', 2995828, 98, 1., 47.13, k_factor=1.),
 
     sample('ZZ',   'ZZ', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 990064, 94, 1.,16.523, k_factor=1.),#NLO from MCFM
 	sample('ZZ_ext',   'ZZ_ext', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM', 998034, 94, 1.,16.523, k_factor=1.),
 
 	sample('WWinclusive', 'WWinclusive', '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 1999000, 208, 1., 12.178, k_factor=1.),#already NNLO xs
 	sample('WW200to600', 'WW200to600', '/WWTo2L2Nu_Mll_200To600_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 1.385, k_factor=1.),#already NNLO xs
     sample('WW600to1200', 'WW600to1200', '/WWTo2L2Nu_Mll_600To1200_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 0.0566, k_factor=1.),#already NNLO xs
     sample('WW1200to2500', 'WW1200to2500', '/WWTo2L2Nu_Mll_1200To2500_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 0.003557, k_factor=1.),#already NNLO xs
     sample('WW2500', 'WW2500','/WWTo2L2Nu_Mll_2500ToInf_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 38969, 208, 1., 0.00005395, k_factor=1.),#already NNLO xs
 
 	sample('Wjets', 'Wjets', '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',29705748,52,1.,61526.7,k_factor=1),#already NNLO xs
 
# 	########sample('ttbar_lep',     'ttbar_lep', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 79092400, 4 , 1., 87.31, k_factor=1.),
 	sample('ttbar_lep50to500',     'ttbar_lep50to500', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 79092400, 4 , 1., 87.31, k_factor=1.),
 	sample('ttbar_lep_500to800',     'ttbar_lep_500to800', '/TTToLL_MLL_500To800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 4 , 1., 0.32611, k_factor=1.),
 	sample('ttbar_lep_800to1200',     'ttbar_lep_800to1200', '/TTToLL_MLL_800To1200_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 199800, 4 , 1., 0.03265, k_factor=1.),
 	sample('ttbar_lep_1200to1800',     'ttbar_lep_1200to1800', '/TTToLL_MLL_1200To1800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 4 , 1., 0.00305, k_factor=1.),
 	sample('ttbar_lep1800toInf',     'ttbar_lep1800toInf', '/TTToLL_MLL_1800ToInf_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 40829, 4 , 1., 0.00017, k_factor=1.),
 	
 	sample('Wantitop', 'WantiTop', '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM',6933094,63 , 1., 35.6, k_factor=1.),#already NNLO xs          
     sample('tW',     'tW', '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM',6952830,66 , 1., 35.6, k_factor=1.),#already NNLO xs
 
# ###    sample('qcd50to80', 'QCD50to80', '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 9954370,43,1.,1,k_factor=1),
 	sample('qcd80to120', 'QCD80to120', '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6986740,43,1.,2762530,k_factor=1),
     sample('qcd120to170', 'QCD120to170', '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6708572,43,1.,471100,k_factor=1),
    	sample('qcd170to300', 'QCD170to300', '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6958708,43,1.,117276,k_factor=1),
     sample('qcd300to470', 'QCD300to470', '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 4150588,43,1.,7823,k_factor=1),
    	sample('qcd470to600', 'QCD470to600', '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3959986,43,1.,648.2,k_factor=1),
     sample('qcd600to800', 'QCD600to800', '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3896412,43,1.,186.9,k_factor=1),
 	sample('qcd800to1000', 'QCD800to1000', '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3992112,43,10,32.293,k_factor=1),
     sample('qcd1000to1400', 'QCD1000to1400', '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',2999069,43,1.,9.4183,k_factor=1),
     sample('qcd1400to1800', 'QCD1400to1800', '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',396409,43,1.,0.84265,k_factor=1),
     sample('qcd1800to2400', 'QCD1800to2400', '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',397660,43,1.,0.114943,k_factor=1),
     sample('qcd2400to3200', 'QCD2400to3200', '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',399226,43,1.,0.00682981,k_factor=1),
     sample('qcd3200', 'QCD3200', '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM',391735,43,1.,0.000165445,k_factor=1),
 
# ### 	N_EVENT scaled by: -N_EVENT * n_neg/n + n * n_pos/n (N_EVENT from report = 28968252; n from weight = 9112991 n_neg = 1507200 (0.1654); n_pos = 7605790 (0.8346); )
 	sample('dyInclusive50', 'DYInclusive50', '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 19385554, 209 , 1., 5765.4,    k_factor=1., is_madgraph=True),  

    ]

samples.reverse()


#if miniAOD:

#else:
for sample in samples:
   exec '%s = sample' % sample.name
   if not miniAOD:
       sample.ana_dataset = '/%s/rradogna-datamc_%s-c4b4ec8fa143ea00cec443e9d0afb38f/USER'  % (sample.dataset.split('/')[1], sample.name)
   else:
       sample.ana_dataset = '/'+ sample.dataset.split('/')[1]+ '/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#dy100to200.ana_dataset = '/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy100to200-1e36332d8badf10b79a5027340f46eb1/USER'
#dy200to400.ana_dataset = '/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy200to400-1e36332d8badf10b79a5027340f46eb1/USER'
#dy400to500.ana_dataset = '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy400to500-1e36332d8badf10b79a5027340f46eb1/USER'
#dy500to700.ana_dataset = '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy500to700-1e36332d8badf10b79a5027340f46eb1/USER'
#dy700to800.ana_dataset = '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy700to800-1e36332d8badf10b79a5027340f46eb1/USER'
#dy800to1000.ana_dataset = '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy800to1000-1e36332d8badf10b79a5027340f46eb1/USER'
#dy1000to1500.ana_dataset = '/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy1000to1500-1e36332d8badf10b79a5027340f46eb1/USER'
#dy1500to2000.ana_dataset = '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy1500to2000-1e36332d8badf10b79a5027340f46eb1/USER'
#dy2000to3000.ana_dataset = '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/rradogna-datamc_dy2000to3000-1e36332d8badf10b79a5027340f46eb1/USER'

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
