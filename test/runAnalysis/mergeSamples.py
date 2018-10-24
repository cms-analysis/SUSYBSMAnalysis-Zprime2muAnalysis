import sys, os, subprocess
sys.path.append('../')
def main():
	from sys import argv
	mc = True
	if argv[1] == "data":
		#samples = ['ana_datamc_Run2017MuonsOnly_SingleMuonRun2017F-ReReco-v1']
		samples = ['ana_datamc_Run2016MuonsOnly_SingleMuonRun2016B-ReReco-v3','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016C-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016D-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016E-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016F-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016G-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v2','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v3','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017B-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017C-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017D-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017E-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017F-ReReco-v1']
		mc = False

	elif argv[1] == "dataElectron":
		#samples = ['ana_datamc_Run2017MuonsOnly_SingleMuonRun2017F-ReReco-v1']
#		samples = ['ana_datamc_Run2016MuonsOnly_SingleMuonRun2016B-ReReco-v3','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016C-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016D-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016E-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016F-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016G-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v2','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v3','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017B-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017C-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017D-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017E-ReReco-v1','ana_datamc_Run2017MuonsOnly_SingleMuonRun2017F-ReReco-v1']
		samples = ['ana_datamc_Run2017_DoubleEG2017B-17Nov2017-v1','ana_datamc_Run2017_DoubleEG2017C-17Nov2017-v1','ana_datamc_Run2017_DoubleEG2017D-17Nov2017-v1','ana_datamc_Run2017_DoubleEG2017E-17Nov2017-v1','ana_datamc_Run2017_DoubleEG2017F-17Nov2017-v1']
		mc = False
	elif argv[1] == "mc2016":
		from histos_resolutionMC import samples2016 as samples
	else:
		from histos_resolutionMC import samples as samples
	if not mc:	
		for sample in samples:
			dirName = sample
		#	dirName = sample.split("/")[1]
			if os.path.isdir(dirName):
				
				fileList  = [dirName+"/"+f for f in os.listdir(dirName) if os.path.isfile(dirName + "/" + f)]
				print "merging %d files for %s"%(len(fileList),sample)
				command = ["hadd","-f","ana_datamc_%s.root"%sample]			
				command += fileList
				subprocess.call(command,stdout=open(os.devnull, 'wb'))

			else:
				print "no output for sample ", dirName

	else:
#		samples = [
 #                      ('dyPt100To250','/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext5-v1/MINIAODSIM'),
#]
		print samples
		for entry in samples:
			dirName = entry[1].split("/")[1]
			sampleName = entry[0]				
			if os.path.isdir(dirName):
				
				fileList  = [dirName+"/"+f for f in os.listdir(dirName) if os.path.isfile(dirName + "/" + f)]
				print "merging %d files for %s"%(len(fileList),sampleName)
				command = ["hadd","-f","ana_datamc_%s.root"%sampleName]			
				command += fileList
				subprocess.call(command,stdout=open(os.devnull, 'wb'))

			else:
				print "no output for sample ", dirName


main()
