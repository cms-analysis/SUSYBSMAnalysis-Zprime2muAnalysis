#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn, files_from_dbs
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import dataset_from_publish_log

miniAOD = True

class sample(object):
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, cross_section_uncert=0.0, k_factor=1, frac_neg_weight=0.0, filenames=None, scheduler='condor', hlt_process_name='HLT', ana_dataset=None, is_madgraph=False, is_zprime=False):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
        self.syst_frac = syst_frac
        self.cross_section = cross_section
        self.cross_section_uncert = cross_section_uncert
        self.frac_neg_weight = frac_neg_weight
        self.k_factor = k_factor
        self.filenames_ = filenames
        self.scheduler = scheduler
        self.hlt_process_name = hlt_process_name
        self.ana_dataset = ana_dataset
        self.is_madgraph = is_madgraph
        self.is_zprime = is_zprime
        # Effective number of events in case of negatively weighted events
        # Slide 32: 
        # https://indico.cern.ch/event/790963/contributions/3306139/attachments/1791391/2919088/CMSWeek201902.pdf
        self.nevents_eff = self.nevents * pow((1-2*self.frac_neg_weight),2)

    # the total weight is partial_weight * integrated_luminosity
    @property
    def partial_weight(self):
        return self.cross_section / float(self.nevents) * self.k_factor 

    @property
    def partial_weight_eff(self):
        return self.cross_section / float(self.nevents_eff) * self.k_factor

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

# XS, cross_section_uncert, frac_neg_weight for all MC samples taken from XSDB which is linked from the DAS page of each dataset
# For example Drell-Yan samples:
# https://cms-gen-dev.cern.ch/xsdb/?columns=57394461&currentPage=0&ordDirection=1&ordFieldName=process_name&pageSize=10&searchQuery=DAS%3DZToMuMu_NNPDF31_13TeV-powheg_M_%2A_%2A
samples = [
    sample('dy50to120', 'DY50to120', '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 2982000, 209 , 1., 2113.0, cross_section_uncert=0.9976,  k_factor=1., frac_neg_weight=0.00974),
    sample('dy120to200', 'DY120to200', '/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 210, 1., 20.55, cross_section_uncert=0.01372, k_factor=1., frac_neg_weight=0.00497),
    sample('dy200to400', 'DY200to400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 211, 1., 2.886, cross_section_uncert=0.001993, k_factor=1., frac_neg_weight=0.00194),
    sample('dy400to800', 'DY400to800', '/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 212, 1., 0.2513, cross_section_uncert=0.000179, k_factor=1., frac_neg_weight=0.00048),
    sample('dy800to1400', 'DY800to1400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', 100000, 72, 1., 0.01707, k_factor=1., cross_section_uncert=1.253e-05, frac_neg_weight=0.00017),
    sample('dy1400to2300', 'DY1400to2300', '/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 70, 1., 0.001366, k_factor=1., cross_section_uncert=1.013E-06, frac_neg_weight=9E-05),
    sample('dy2300to3500', 'DY2300to3500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 70, 1.,8.178e-05, k_factor=1., cross_section_uncert=5.997e-08, frac_neg_weight=0.00015),
    sample('dy3500to4500', 'DY3500to4500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 70, 1., 3.191e-06, k_factor=1., cross_section_uncert=2.243e-09, frac_neg_weight=0.00101),
    sample('dy4500to6000', 'DY4500to6000', '/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 70, 1., 2.787e-07, k_factor=1., cross_section_uncert=1.702e-10, frac_neg_weight=0.00349),
    sample('dy6000toInf', 'DY6000toInf', '/ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 100000, 70, 1., 9.569e-09, k_factor=1., cross_section_uncert=8.431e-12, frac_neg_weight=0.00933),

    # Diboson
    sample('WZ', 'WZ', '/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM', 3885000, 98, 1., 27.6, k_factor=1., cross_section_uncert=0.04),
    #sample('WZ', 'WZ', '/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', XXXX, 98, 1., 27.6, k_factor=1., cross_section_uncert=0.04), # in production
    sample('ZZ', 'ZZ', '/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 1979000, 94, 1., 12.14, k_factor=1., cross_section_uncert=0.01964),

    sample('WW', 'WW', '/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', 3885000, 208, 1., 75.8, k_factor=1., cross_section_uncert=0.1123),
    #sample('WW', 'WW', '/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', XXXX, 208, 1., 75.8, k_factor=1., cross_section_uncert=0.1123), # in production

    # WW Samples still in production
    #sample('WWinclusive', 'WWinclusive', '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 1999000, 208, 1., 12.178, k_factor=1.),#already NNLO xs
    #sample('WW200to600', 'WW200to600', '/WWTo2L2Nu_Mll_200To600_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 1.385, k_factor=1.),#already NNLO xs
    #sample('WW600to1200', 'WW600to1200', '/WWTo2L2Nu_Mll_600To1200_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 0.0566, k_factor=1.),#already NNLO xs
    #sample('WW1200to2500', 'WW1200to2500', '/WWTo2L2Nu_Mll_1200To2500_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 208, 1., 0.003557, k_factor=1.),#already NNLO xs
    #sample('WW2500', 'WW2500','/WWTo2L2Nu_Mll_2500ToInf_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 38969, 208, 1., 0.00005395, k_factor=1.),#already NNLO xs
 
 	sample('Wjets', 'Wjets', '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM', 60750603, 52, 1., 52940.0, k_factor=1, cross_section_uncert=121.9, frac_neg_weight=0.0004079), # Usable for now? Switch to recommended Autumn18MiniAOD campaign
    #sample('Wjets', 'Wjets', '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM', XXXX, 52, 1., 52940.0, k_factor=1, cross_section_uncert=121.9, frac_neg_weight=0.0004079), # This one is the recommended MC Campaign
 
    # t-tbar samples binned in dilepton mass still in production
    sample('ttbar_lep', 'ttbar', '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', 64310000, 4 , 1., 88.29, k_factor=1.),
    #sample('ttbar_lep50to500',     'ttbar_lep50to500', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 79092400, 4 , 1., 87.31, k_factor=1.),
    #sample('ttbar_lep_500to800',     'ttbar_lep_500to800', '/TTToLL_MLL_500To800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 4 , 1., 0.32611, k_factor=1.),
 	#sample('ttbar_lep_800to1200',     'ttbar_lep_800to1200', '/TTToLL_MLL_800To1200_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 199800, 4 , 1., 0.03265, k_factor=1.),
 	#sample('ttbar_lep_1200to1800',     'ttbar_lep_1200to1800', '/TTToLL_MLL_1200To1800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 200000, 4 , 1., 0.00305, k_factor=1.),
 	#sample('ttbar_lep1800toInf',     'ttbar_lep1800toInf', '/TTToLL_MLL_1800ToInf_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 40829, 4 , 1., 0.00017, k_factor=1.),
 	
    # Single-Top
 	sample('tbarW', 'tbarW', '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM', 8000000, 63 , 1., 34.97, k_factor=1., cross_section_uncert=0.02827, frac_neg_weight=0.0034),
 	sample('tbarW_ext', 'tbarW', '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', 7623000, 63 , 1., 34.97, k_factor=1., cross_section_uncert=0.02827, frac_neg_weight=0.0034),
     sample('tW', 'tW', '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM', 9438000, 66 , 1., 34.91, k_factor=1., cross_section_uncert=0.02817, frac_neg_weight=0.003758),
     sample('tW_ext', 'tW', '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM', 6952830, 66 , 1., 34.91, k_factor=1.,cross_section_uncert=0.02817, frac_neg_weight=0.003758),
     
    # QCD Samples
    #sample('qcd50to80', 'QCD50to80', '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 9954370,43,1.,1,k_factor=1),
    #sample('qcd80to120', 'QCD80to120', '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6986740,43,1.,2762530,k_factor=1),
    #sample('qcd120to170', 'QCD120to170', '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6708572,43,1.,471100,k_factor=1),
    #sample('qcd170to300', 'QCD170to300', '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6958708,43,1.,117276,k_factor=1),
    #sample('qcd300to470', 'QCD300to470', '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 4150588,43,1.,7823,k_factor=1),
    #sample('qcd470to600', 'QCD470to600', '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3959986,43,1.,648.2,k_factor=1),
    #sample('qcd600to800', 'QCD600to800', '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3896412,43,1.,186.9,k_factor=1),
    #sample('qcd800to1000', 'QCD800to1000', '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3992112,43,10,32.293,k_factor=1),
    #sample('qcd1000to1400', 'QCD1000to1400', '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',2999069,43,1.,9.4183,k_factor=1),
    #sample('qcd1400to1800', 'QCD1400to1800', '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',396409,43,1.,0.84265,k_factor=1),
    #sample('qcd1800to2400', 'QCD1800to2400', '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',397660,43,1.,0.114943,k_factor=1),
    #sample('qcd2400to3200', 'QCD2400to3200', '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',399226,43,1.,0.00682981,k_factor=1),
    #sample('qcd3200', 'QCD3200', '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM',391735,43,1.,0.000165445,k_factor=1),
 
# ### 	N_EVENT scaled by: -N_EVENT * n_neg/n + n * n_pos/n (N_EVENT from report = 28968252; n from weight = 9112991 n_neg = 1507200 (0.1654); n_pos = 7605790 (0.8346); )

    # DYJetsToLL - in the past we've used amcatnloFXFX but why?
    #sample('dyInclusive50_amcatnlo', 'DYInclusive50', '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', 997561, 209 , 1., 6529.0, k_factor=1., is_madgraph=True, cross_section_uncert=28.29, frac_neg_weight=0.1624),  
 	#sample('dyInclusive50_madgraph', 'DYInclusive50', '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM', 100194597, 209 , 1., 5343.0, k_factor=1., is_madgraph=True, cross_section_uncert=12.64, frac_neg_weight=0.0004962), # NNLO XS is 6225.42 and has fractional uncertainty of 2% see 
    # https://indico.cern.ch/event/746829/contributions/3138541/attachments/1717905/2772129/Drell-Yan_jets_crosssection.pdf

    ]

samples.reverse()

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
