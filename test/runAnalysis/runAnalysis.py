import argparse, subprocess, os	



def getFilterSnippet(name):
	print name
	ttbarFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
			       src = cms.InputTag('prunedGenParticles'),
			       min_mass = cms.double(50),
			       max_mass = cms.double(500), 
			       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	wwFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
				       src = cms.InputTag('prunedGenParticles'),
				       min_mass = cms.double(50),
				       max_mass = cms.double(200), 
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	     
	dyFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
				       src = cms.InputTag('prunedGenParticles'),                                      
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	if "dyInclusive" in name:
		return dyFilter
#	elif "ttbar_lep50to500" in name:
#		return ttbarFilter
#	elif "WWinclusive" in name:
#		return wwFilter
	else:
		return ""

def getCRABCfgWeirdSubmission(name,dataset,fileList):

	crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_%s'
config.General.workArea = 'crab'
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'
#config.JobType.priority = 1
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_%s'
config.Data.outLFNDirBase = '/store/user/jschulte'
config.Site.storageSite = 'T2_US_Purdue'
config.Site.whitelist = ['T2_ES_IFCA','T2_US_Nebraska','T2_US_UCSD']
#config.Site.whitelist = ['T1_US_FNAL']
config.Data.userInputFiles = %s
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB  = 8000
'''
	result = crab_cfg%(name,dataset,name,fileList)
	
	return result




def getCRABCfg(name,dataset,lumi_mask=""):

	crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_%s'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '%s'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_%s'
config.Data.outLFNDirBase = '/store/user/jschulte'
#config.Data.ignoreLocality = True 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
config.JobType.maxMemoryMB  = 8000
%s
'''
	data_config='''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
config.Data.lumiMask = '%s'
'''

	mc_config='''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 400000
'''
	if lumi_mask =="":
		result = crab_cfg%(name,dataset,name,mc_config)
	else:
		data_cfg = data_config%lumi_mask	
		result = crab_cfg%(name,dataset,name,data_cfg)
	
	return result




def main():

	parser = argparse.ArgumentParser(description='tool to run dilepton analysis')
	
	parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
						  help="Verbose mode.")
	parser.add_argument("-d", "--data", action="store_true", dest="data", default=False,help="run on data")
	parser.add_argument("-l", "--local", action="store_true", dest="local", default=False,help="run locally")
	parser.add_argument("-s", "--submit", action="store_true", dest="submit", default=False,help="submit to CRAB")
	parser.add_argument("-r", "--resolution", action="store_true", dest="resolution", default=False,help="run jobs for resolution studies")
	parser.add_argument("-w", "--write", action="store_true", dest="write", default=False,help="write config but not execute")
	parser.add_argument("-e", "--electrons", action="store_true", dest="electrons", default=False,help="run electrons")
	parser.add_argument( "--ci2017", action="store_true", dest="ci2017", default=False,help="run CI MC for 2017")
	parser.add_argument( "--2016", action="store_true", dest="do2016", default=False,help="run for 2016")
	parser.add_argument( "--2018", action="store_true", dest="do2018", default=False,help="run for 2018")
	args = parser.parse_args()
	isMC = "True"
	GT = "94X_mc2017_realistic_v14"
	if args.data:
		GT = "94X_dataRun2_ReReco_EOY17_v6"
		isMC = 'False'
	arguments = {}
	arguments["GT"] = GT
	arguments["isMC"] = isMC
	cmssw_cfg = open('setup.py').read()%arguments
	prefix = "muons_"	
	if not args.resolution:
		
		if args.electrons:
			prefix = "electrons_"
			cmssw_cfg += open('histogrammerElectrons.py').read()
		else:	
			cmssw_cfg += open('histogrammerMuons.py').read()
	else:
		prefix = "resolution_"
		cmssw_cfg += open('resolution.py').read()
	
	open('cmssw_cfg.py', 'wt').write(cmssw_cfg)

	if not args.write:
		if args.local:
			subprocess.call(['cmsRun','cmssw_cfg.py'])	
		else:
			if args.electrons and args.data:
				if args.do2018:
					from samples import data_electrons_2018 as samples
				elif args.do2016:	
					from samples import data_electrons_2016 as samples
				else:
					from samples import data_electrons_2017 as samples
			elif args.electrons and not args.data:
				if args.ci2017:	
					from samples import ci_electrons_2017 as samples
				else:
					from samples import backgrounds_electrons_2017 as samples
			elif args.data:
				if args.do2018:
					from samples import data_muons_2018 as samples
				elif args.do2016:
					from samples import data_muons_2016 as samples
				else:
					from samples import data_muons_2017 as samples
			else:
				if args.ci2017:
					from samples import ci_muons_2017 as samples 
				else:	
					from samples import backgrounds_muons_2017 as samples 
			lumi_mask = ""
			GT = "94X_mc2017_realistic_v14"
			if args.data:
				if args.electrons: 
					lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
					if args.do2018:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-324420_13TeV_PromptReco_Collisions18_JSON.txt"
					if args.do2016:	
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

				else:
					lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt"
					if args.do2018:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-324420_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt"
					if args.do2016:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_MuonPhys.txt"
				GT = "94X_dataRun2_ReReco_EOY17_v6"
			for dataset_name,  dataset in samples:
				cmssw_tmp = cmssw_cfg
			 	crab_cfg = getCRABCfg(prefix+dataset_name,dataset,lumi_mask)
		                open('crabConfig.py', 'wt').write(crab_cfg)
				cmssw_tmp+=getFilterSnippet(dataset_name)
				if "dy" in dataset_name:
					if "HistosFromPAT.usekFactor = False" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = False','HistosFromPAT.usekFactor = True')
				else:
					if "HistosFromPAT.usekFactor = True" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = True','HistosFromPAT.usekFactor = False')
				#print getFilterSnippet(dataset_name)
				open('cmssw_cfg.py', 'wt').write(cmssw_tmp)
            			if args.submit:
                			os.system('crab submit -c crabConfig.py')

			if args.resolution and not args.data:
				print "submitting also weird samples"
				from samples import resolution_extra as samples2
			
				from dbs.apis.dbsClient import DbsApi
				dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
			       
				for name, ana_dataset in samples2:
				    cmssw_tmp = cmssw_cfg
				    fileDictList=dbs.listFiles(dataset=ana_dataset)
				    
				    print ("dataset %s has %d files" % (ana_dataset, len(fileDictList)))
				   
				# DBS client returns a list of dictionaries, but we want a list of Logical File Names
				    lfnList = [ dic['logical_file_name'] for dic in fileDictList ]	
				    crab_cfg = getCRABCfgWeirdSubmission(prefix+dataset_name,dataset,lfnList)
				    open('crabConfig.py', 'wt').write(crab_cfg)
				    #cmssw_cfg+=getFilterSnippet(dataset_name) # high mass tails not available yet
				    open('cmssw_cfg.py', 'wt').write(cmssw_tmp)
				    if args.submit:
					os.system('crab submit -c crabConfig.py')


main()	
