import argparse, subprocess, os, sys



def getFilterSnippet(name,args,doApply=False,year=2017):
	print name
	ttbarFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('TTbarSelection',
			       src = cms.InputTag('prunedGenParticles'),
			       min_mass = cms.double(50),
			       max_mass = cms.double(500), 
			       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''

	wwFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
				       src = cms.InputTag('prunedGenParticles'),
				       min_mass = cms.double(50),
				       max_mass = cms.double(200), 
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''

	     
	dyFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
				       src = cms.InputTag('prunedGenParticles'),                                      
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''
	CIFilter300 = '''
process.DYGenMassFilter = cms.EDFilter('QScaleSelector',
				       src = cms.InputTag('generator'),
				       min_mass = cms.double(300),
				       max_mass = cms.double(800),
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''
	CIFilter800 = '''
process.DYGenMassFilter = cms.EDFilter('QScaleSelector',
				       src = cms.InputTag('generator'),
				       min_mass = cms.double(800),
				       max_mass = cms.double(1300),
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''
	CIFilter1300 = '''
process.DYGenMassFilter = cms.EDFilter('QScaleSelector',
				       src = cms.InputTag('generator'),
				       min_mass = cms.double(1300),
				       max_mass = cms.double(2000),
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''
	ZPtFilter = '''    
process.DYGenMassFilter = cms.EDFilter('DyPt_ZSkim',
                                       src = cms.InputTag('prunedGenParticles'),
                                       min_mass = cms.double(0),
                                       max_mass = cms.double(50), 
                                       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(2,process.DYGenMassFilter)'''

	if "dyInclusive" in name:
		if args.resolution:
			return ZPtFilter 
		else:
			return dyFilter
	elif (name == "ttbar_lep" or name == "ttbar_lep_ext" or name == "ttbar_lep50to500") and doApply:
		return ttbarFilter
	elif "WWinclusive" in name and doApply:
		return wwFilter
	elif "CI" in name and year==2016:
		if "M300" in name:
			return CIFilter300
		elif "M800" in name:
			return CIFilter800
		elif "M1300" in name and "ConLL" in name:
			return CIFilter1300
		else:
			return ""
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
#config.Site.whitelist = ['T2_ES_IFCA','T2_US_MIT','T2_US_UCSD']
config.Data.userInputFiles = %s
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.JobType.maxMemoryMB  = 8000
config.JobType.allowUndistributedCMSSW = True
'''
	result = crab_cfg%(name,name,fileList)
	
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
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = 'dileptonAna_%s'
config.Data.outLFNDirBase = '/store/user/jschulte/'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True
config.Site.blacklist = ['T2_US_Caltech']
%s
'''
	data_config='''
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.lumiMask = '%s'
'''

	mc_config='''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000
'''
	if lumi_mask =="":
		result = crab_cfg%(name,dataset,name,mc_config)
	else:
		data_cfg = data_config%lumi_mask	
		result = crab_cfg%(name,dataset,name,data_cfg)
	
	return result


def getCRABCfgAAA(name,dataset,lumi_mask=""):

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
config.Data.outLFNDirBase = '/store/user/jschulte/'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
#config.General.instance = 'preprod' 
config.Site.whitelist = ["T2_US_*"]
config.Site.blacklist = ['T2_US_Caltech']
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True
%s
'''
	data_config='''
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.lumiMask = '%s'
'''

	mc_config='''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000
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
	#parser.add_argument("-l", "--local", action="store_true", dest="local", default=False,help="run locally")
	#parser.add_argument("--allLocal", action="store_true", dest="allLocal", default=False,help="run a full list of samples locally")
	parser.add_argument("--file", action="store", dest="inputFile", type=str  , default="",help="file to run over locally")
	parser.add_argument("--dataset", action="store", dest="inputDataset", type=str , default="",help="dataset to run over locally")
	parser.add_argument("-s", "--submit", action="store_true", dest="submit", default=False,help="submit to CRAB")
	parser.add_argument("-r", "--resolution", action="store_true", dest="resolution", default=False,help="run jobs for resolution studies")
	parser.add_argument("--pdf", action="store_true", dest="pdf", default=False,help="produce trees for PDF studies")
	parser.add_argument("--genMassOther", action="store_true", dest="genMassOther", default=False,help="produce gen mass histograms for other backgrounds")
	parser.add_argument("--genMass", action="store_true", dest="genMass", default=False,help="produce trees for genMas studies")
	parser.add_argument("-w", "--write", action="store_true", dest="write", default=False,help="write config but not execute")
	parser.add_argument("-e", "--electrons", action="store_true", dest="electrons", default=False,help="run electrons")
	parser.add_argument( "--ci2016", action="store_true", dest="ci2016", default=False,help="run CI MC for 2016")
	parser.add_argument( "--ci2017", action="store_true", dest="ci2017", default=False,help="run CI MC for 2017")
	parser.add_argument( "--ci2018", action="store_true", dest="ci2018", default=False,help="run CI MC for 2018")
	parser.add_argument( "--2016", action="store_true", dest="do2016", default=False,help="run for 2016")
	parser.add_argument( "--2018", action="store_true", dest="do2018", default=False,help="run for 2018")
	parser.add_argument( "--addNTuples", action="store_true", dest="addNTuples", default=False,help="add nTuples to histogrammer workflow")
	parser.add_argument( "--add2016", action="store_true", dest="add2016", default=False, help="run ADD MC for 2016")
	parser.add_argument( "--add2017", action="store_true", dest="add2017", default=False, help="run ADD MC for 2017")
	parser.add_argument( "--add2018", action="store_true", dest="add2018", default=False, help="run ADD MC for 2018")
	args = parser.parse_args()



	if args.resolution and args.addNTuples:
		print "warning, addNTuplets does nothing for resolution workflow"

	if args.ci2018 or args.add2018:
		args.do2018 = True

	if args.ci2016 or args.add2016:
		args.do2016 = True

	isMC = "True"
	#GT = "94X_mc2017_realistic_v14"
	GT = "94X_mc2017_realistic_v17"
	if args.data:
		GT = "94X_dataRun2_v11"
		isMC = 'False'
	arguments = {}
	arguments["GT"] = GT
	arguments["isMC"] = isMC
	arguments["addNTuples"] = args.addNTuples
	arguments["year"] = 2017
	if args.do2016:
		arguments["year"] = 2016
		GT = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
		if not isMC:
			GT = '80X_dataRun2_2016SeptRepro_v6'
	if args.do2018:
		arguments["year"] = 2018
		GT = '102X_upgrade2018_realistic_v18'
		if not isMC:
			GT = '102X_dataRun2_Sep2018ABC_v2'
	
	prefix = "muons_"	
	cmssw_cfg = open('setup.py').read()

	if args.resolution:
		prefix = "resolution_"
		analyzer= open('resolution.py').read()
	elif args.pdf:
		prefix = "pdf_"
		analyzer= open('pdfStudies.py').read()

	elif args.genMass:
		prefix = "genMass_"
		if args.electrons:
			prefix += "electrons_"
		analyzer= open('genMass.py').read()
	elif args.genMassOther:
		prefix = "genMassOther_"
		if args.electrons:
			prefix += "electrons_"	
		analyzer= open('genMassOther.py').read()

	else:
		
		if args.electrons:
			prefix = "electrons_"
			analyzer= open('histogrammerElectrons.py').read()
		else:	
			analyzer= open('histogrammerMuons.py').read()
	applyAllGenFilters = True
	if args.do2016:
		prefix = prefix + "2016_"
	if args.do2018:
		prefix = prefix + "2018_"	
#	open('cmssw_cfg.py', 'wt').write(cmssw_cfg)
#	print prefix
	if not args.write:
		#if args.local:
		#	subprocess.call(['cmsRun','cmssw_cfg.py'])	
		#else:
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
				elif args.ci2018:	
					from samples import ci_electrons_2017 as samples
				elif args.ci2016:	
					from samples import ci_electrons_2016 as samples
				elif args.add2017:	
					from samples import add_2017 as samples
				elif args.add2018:	
					from samples import add_2018 as samples
				elif args.add2016:	
					from samples import add_2016 as samples
				elif args.do2016:
					from samples import backgrounds_electrons_2016 as samples			
				elif args.do2018:
					from samples import backgrounds_electrons_2018 as samples
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
				elif args.ci2018:
					from samples import ci_muons_2017 as samples 
				elif args.ci2016:
					from samples import ci_muons_2016 as samples 
				elif args.add2018:
					from samples import add_2017 as samples 
				elif args.add2017:
					from samples import add_2017 as samples 
				elif args.add2016:
					from samples import add_2016 as samples 
				elif args.resolution:
					if args.do2016:
						from samples import resolution_2016 as samples 
					elif args.do2018:
						from samples import resolution_2018 as samples 
					else:	
						from samples import resolution_2017 as samples 
				elif args.do2016:	
					from samples import backgrounds_muons_2016 as samples			
				elif args.do2018:	
					from samples import backgrounds_muons_2018 as samples
				else:	
					from samples import backgrounds_muons_2017 as samples 

			if not args.inputFile == "":
				print "running over single file %s"%args.inputFile
				samples = [("dummy", args.inputFile)]
				if args.submit:
					print "can't submit jobs over a single file, switching to local"
					args.submit = False
			elif not args.inputDataset == "":
				validSample = False
				for sample in samples:
					if args.inputDataset == sample[1]:
						validSample = True
						samples = [sample]
				if validSample:
					print "running over sample %s"%args.inputDataset
				else:
					print "sample does not fit the chosen configuration, please reconsider"
					sys.exit()


			lumi_mask = ""
			GT = "94X_mc2017_realistic_v14"
			if args.add2016:
				GT = "80X_mcRun2_asymptotic_2016_TrancheIV_v6"
			if args.data:
				if args.electrons: 
					lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
					if args.do2018:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
					if args.do2016:	
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

				else:
					lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt"
					if args.do2018:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt"
					if args.do2016:
						lumi_mask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_MuonPhys.txt"
				GT = "94X_dataRun2_ReReco_EOY17_v6"


#			lumi_mask = '/afs/cern.ch/work/j/jschulte/test/CMSSW_10_2_15_patch1/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/runAnalysis/crab/crab_dileptonAna_muons_2016_SingleMuonRun2016B-23Sep2016_v3/results/notFinishedLumis.json'
			for dataset_name,  dataset in samples:
				#print dataset_name, "bkub"
				#//if not ("LR" in dataset_name or "RL" in dataset_name): continue	
				if args.do2018 and args.electrons and dataset_name == "dy50to120":
					lumi_mask="dy2018JSON.txt"
				
				arguments['name'] = dataset_name
				cmssw_tmp = cmssw_cfg
				cmssw_tmp = cmssw_tmp%arguments
				cmssw_tmp += analyzer
				
				cmssw_tmp+=getFilterSnippet(dataset_name,args,applyAllGenFilters,year=arguments["year"])
				if args.genMassOther:
					cmssw_tmp = cmssw_tmp.replace('getattr(process,path_name).insert(2,process.DYGenMassFilter)','getattr(process,path_name).insert(1,process.DYGenMassFilter)')
				if "dy" in dataset_name:
					if "HistosFromPAT.usekFactor = False" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = False','HistosFromPAT.usekFactor = True')
				else:
					if "HistosFromPAT.usekFactor = True" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = True','HistosFromPAT.usekFactor = False')
				if "ttbar" in dataset_name and not args.do2016:
					if "HistosFromPAT.useTTBarWeight = False" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.useTTBarWeight = False','HistosFromPAT.useTTBarWeight = True')
				else:
					if "HistosFromPAT.useTTBarWeight = True" in cmssw_tmp:
						cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.useTTBarWeight = True','HistosFromPAT.useTTBarWeight = False')

				if not args.do2016 and not args.ci2017 and not args.ci2016 and not args.ci2018 and not args.add2016 and not args.add2017 and not args.add2018 and not args.data and not args.resolution and args.electrons and not dataset_name == 'dummy':
					print "trying"
					if args.do2018:
						cmssw_tmp = cmssw_tmp.replace('mc_2018',dataset_name)
					else:
						cmssw_tmp = cmssw_tmp.replace('mc_2017',dataset_name)

				if args.ci2018 or args.add2018:
						cmssw_tmp = cmssw_tmp.replace('mc_2018','mc_2017')
				#if args.do2018 and 'dy' in dataset_name and not dataset_name == "dyInclusive50":
				#	cmssw_tmp = cmssw_tmp.replace('mc_2018','mc_2018_flat')
				#toReweight = ['WW200to600','WW600to1200','WW1200to2500','WW2500',"ttbar_lep_ext","ttbar_lep",'ttbar_lep_500to800_ext','ttbar_lep_500to800','ttbar_lep_800to1200','ttbar_lep_1200to1800','ttbar_lep1800toInf']
				#if (args.do2018 and dataset_name in toReweight) or (args.do2018 and 'CITo2E' in dataset_name):
				#	cmssw_tmp = cmssw_tmp.replace('mc_2018','mc_2017')

				if args.submit:
					os.system('dasgoclient --query="site dataset=%s" > sites.txt'%dataset)

					with open("sites.txt") as f:
						content = f.readlines()
					content = [x.strip() for x in content] 
					nT2 = 0
					for site in content:
						if "T2" in site:
							nT2 += 1				
					if nT2 <= 1:
						crab_cfg = getCRABCfgAAA(prefix+dataset_name,dataset,lumi_mask)
					else:	
						crab_cfg = getCRABCfg(prefix+dataset_name,dataset,lumi_mask)
					if args.do2016 and ("Lam100kTeV" in dataset_name or "Lam10TeV"	in dataset_name):
						crab_cfg = crab_cfg.replace("config.Data.inputDBS = 'global'","config.Data.inputDBS = 'phys03'")
					if args.do2016 and dataset_name == "CITo2E_Lam1TeVConLR_M300":
						crab_cfg = crab_cfg + '\n'					
						crab_cfg = crab_cfg + 'config.Data.allowNonValidInputDataset = True'					
					open('crabConfig.py', 'wt').write(crab_cfg)

					#print getFilterSnippet(dataset_name)
					open('cmssw_cfg.py', 'wt').write(cmssw_tmp)
					os.system('crab submit -c crabConfig.py')


				else:
					# write add filename information
					if dataset_name == "dummy":
						cmssw_tmp = cmssw_tmp.replace('dummyFile', dataset)
					else:
						if 'CITo2Mu_Lam10TeV' in dataset_name and args.ci2016:
							print 'here'	
							os.system('dasgoclient -query="file dataset=%s instance=prod/phys03| grep file.name" > myfiles.txt'%dataset)
						else:	
							os.system('dasgoclient -query="file dataset=%s | grep file.name" > myfiles.txt'%dataset)
						with open('myfiles.txt') as fin:
							content = fin.readlines()
						content = [x.strip() for x in content]
						files = "["
						for fileName in content:
							files += "'%s'"%fileName
							files += ","
						files += "]"
						fin.close()
						cmssw_tmp = cmssw_tmp.replace("['dummyFile']", files)

					open('cmssw_cfg.py', 'wt').write(cmssw_tmp)

					subprocess.call(['cmsRun','cmssw_cfg.py'])	

			if args.do2016 and not args.resolution and not args.data and not args.ci2016 and not args.add2016:
				print "submitting also weird samples"
				from samples import backgrounds_electrons_2016_extra as samples2
			
				from dbs.apis.dbsClient import DbsApi
				dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
			       
				for name, ana_dataset in samples2:
					arguments['name'] = dataset_name
					cmssw_tmp = cmssw_cfg
					cmssw_tmp = cmssw_tmp%arguments
					cmssw_tmp += analyzer
					
					cmssw_tmp+=getFilterSnippet(dataset_name,args,applyAllGenFilters,year=arguments["year"])
					if "dy" in dataset_name:
						if "HistosFromPAT.usekFactor = False" in cmssw_tmp:
							cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = False','HistosFromPAT.usekFactor = True')
					else:
						if "HistosFromPAT.usekFactor = True" in cmssw_tmp:
							cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.usekFactor = True','HistosFromPAT.usekFactor = False')
					if "ttbar" in dataset_name:
						if "HistosFromPAT.useTTBarWeight = False" in cmssw_tmp:
							cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.useTTBarWeight = False','HistosFromPAT.useTTBarWeight = True')
					else:
						if "HistosFromPAT.useTTBarWeight = True" in cmssw_tmp:
							cmssw_tmp = cmssw_tmp.replace('HistosFromPAT.useTTBarWeight = True','HistosFromPAT.useTTBarWeight = False')

					if not args.do2016 and not args.do2018 and not args.ci2017 and not args.data and not args.resolution and args.electrons and not dataset_name == 'dummy':
						print "trying"
						cmssw_tmp = cmssw_tmp.replace('mc_2017',dataset_name)
					#if args.do2018 and 'dy' in dataset_name and not dataset_name == "dyInclusive50":
					#	cmssw_tmp = cmssw_tmp.replace('mc_2018','mc_2018_flat')
					toReweight = ['WW200to600','WW600to1200','WW1200to2500','WW2500',"ttbar_lep_ext","ttbar_lep",'ttbar_lep_500to800_ext','ttbar_lep_500to800','ttbar_lep_800to1200','ttbar_lep_1200to1800','ttbar_lep1800toInf']
					if (args.do2018 and dataset_name in toReweight) or (args.do2018 and 'CITo2E' in dataset_name):
						cmssw_tmp = cmssw_tmp.replace('mc_2018','mc_2017')


				        fileDictList=dbs.listFiles(dataset=ana_dataset)
				    
				        print ("dataset %s has %d files" % (ana_dataset, len(fileDictList)))
				   
				        # DBS client returns a list of dictionaries, but we want a list of Logical File Names
				        lfnList = [ dic['logical_file_name'] for dic in fileDictList ]	
				        crab_cfg = getCRABCfgWeirdSubmission(prefix+name,ana_dataset,lfnList)
				        open('crabConfig.py', 'wt').write(crab_cfg)
				        #cmssw_cfg+=getFilterSnippet(dataset_name) # high mass tails not available yet
				        #cmssw_tmp = cmssw_tmp.replace('/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root', lfnList[0])
				        open('cmssw_cfg.py', 'wt').write(cmssw_tmp)
				        if args.submit:
						os.system('crab submit -c crabConfig.py')
						#print "Test!"



			if args.resolution and not args.data and not args.do2016 and not args.do2018:
				print "submitting also weird samples"
				from samples import resolution_extra as samples2
			
				from dbs.apis.dbsClient import DbsApi
				dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
			       
				for name, ana_dataset in samples2:
				    arguments['name'] = name
				    cmssw_tmp = cmssw_cfg
				    cmssw_tmp = cmssw_tmp%arguments
				    cmssw_tmp += analyzer

				    fileDictList=dbs.listFiles(dataset=ana_dataset)
				    
				    print ("dataset %s has %d files" % (ana_dataset, len(fileDictList)))
				   
				# DBS client returns a list of dictionaries, but we want a list of Logical File Names
				    lfnList = [ dic['logical_file_name'] for dic in fileDictList ]	
				    crab_cfg = getCRABCfgWeirdSubmission(prefix+name,ana_dataset,lfnList)
				    open('crabConfig.py', 'wt').write(crab_cfg)
				    #cmssw_cfg+=getFilterSnippet(dataset_name) # high mass tails not available yet
				    #cmssw_tmp = cmssw_tmp.replace('/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/50000/00105BAD-63E0-E711-8640-02163E0146C5.root', lfnList[0])
				    open('cmssw_cfg.py', 'wt').write(cmssw_tmp)
				    if args.submit:
					os.system('crab submit -c crabConfig.py')
					#print "Test!"


main()
