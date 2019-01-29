from CRABAPI.RawCommand import crabCommand
import subprocess
from glob import glob
import os
# If you want crabCommand to be quiet:
from CRABClient.UserUtilities import setConsoleLogLevel
from CRABClient.ClientUtilities import LOGLEVEL_MUTE
setConsoleLogLevel(LOGLEVEL_MUTE)

import argparse


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

parser = argparse.ArgumentParser(description='tool to run dilepton analysis')
	
parser.add_argument("--forceMerge", action="store_true", dest="forceMerge", default=False, help="force merge even if not all files present")

args = parser.parse_args()

# With this function you can change the console log level at any time.
for d in glob("crab2/*/"):
#	if not "CITo2Mu" in d: 
#^		continue
	print "check status of task %s"%d
	res = crabCommand('status',dir=d)
	if res.get("jobList"):
	    nJobs = len(res.get("jobList"))
	    nFinishedJobs = 0
	    failedJobs = []
	    finishedJobs = []
	    for job in res.get("jobList"):
		if job[0] == 'finished':
			finishedJobs.append(job[1])
		elif job[0] == "failed":
			failedJobs.append(job[1])
	    nFinishedJobs = len(finishedJobs)
	    if nFinishedJobs > 0:
		print "some jobs done, checking output"
		fileCount = len(glob(d+'/results/*.root'))
		if fileCount == nJobs:
			"all files already downloaded, doing nothing"
			continue
		if fileCount == 0:
			"no files downloaded yet, starting"
			if nJobs < 500:
				#subprocess.call(["crab",'getoutput',d])
				resGet = crabCommand("getoutput",dir=d)
			elif nJobs < 1000:
				resGet = crabCommand("getoutput","--jobids", "1-500",dir=d)
				resGet2 = crabCommand("getoutput","--jobids", "501-%d"%nJobs,dir=d)
				#subprocess.call(["crab",'getoutput',d, "--jobids", "1-500"])
				#subprocess.call(["crab",'getoutput',d, "--jobids", "501-%d"%nJobs])i
				
				resGet['failed'].update( resGet2['failed'])
				resGet['success'].update(resGet2['success'])
			else:
				subprocess.call(["crab",'getoutput',d, "--jobids", "1-500"])
				subprocess.call(["crab",'getoutput',d, "--jobids", "501-1000"])
				subprocess.call(["crab",'getoutput',d, "--jobids", "1001-%d"%nJobs])

		else:
			jobIds = ""
			jobsToGet = []
 			for i in range(1,nJobs+1):
				if not os.path.isfile(d+'/results/zp2mu_histos_%d.root'%i) and str(i) in finishedJobs:
					jobsToGet.append(i)
			if len(jobsToGet) == 0:
				print "no new files to get"
				resGet = {}
				resGet["failed"] = {}
				resGet["success"] = {}
			else:
				print jobsToGet
				print "getting %d additional files"%len(jobsToGet)
				if len(jobsToGet) > 500:
					sublists = chunks(jobsToGet,499)
					resGet = {}
					resGet['success'] = []
					resGet['failed'] = []
					for sublist in sublists:
						jobIds = ""
					
						for i in sublist:	
							if jobIds == "":
								jobIds += "%d"%i
							else:
								jobIds += ",%d"%i

						resGetTemp = crabCommand("getoutput","--jobids", jobIds,dir=d)
						resGet['failed'].update( resGetTemp['failed'])
						resGet['success'].update(resGetTemp['success'])
				else:
					for i in jobsToGet:	
						if jobIds == "":
							jobIds += "%s"%i
						else:
							jobIds += ",%d"%i
					resGet = crabCommand("getoutput","--jobids", jobIds,dir=d)
		print "got %d new output files, merging if necessary"%len(resGet['success'])
		if (args.forceMerge or len(finishedJobs) == nJobs ) and (len(resGet['failed']) == 0 and not len(resGet['success']) == 0):
			if nJobs == 1:
				cpCommand = ["cp", d+'/results/zp2mu_histos_1.root', d.split('/')[-2].split('crab_')[-1]+".root"]
				print "only one output file, copying instead of merging"
				subprocess.call(cpCommand)
			else:
				haddCommand = ['hadd',"-f",d.split('/')[-2].split('crab_')[-1]+".root"]
		
				for i in finishedJobs:
					haddCommand.append(d+'/results/zp2mu_histos_%s.root'%i)		
				subprocess.call(haddCommand)
		else:
			print "some downloads failed or no new files found - not merging the output"
	    if len(failedJobs) != 0:
		print "some jobs failed, resubmitting"
		jobIds = ""
		nToResubmit = 0
		resResubmit = crabCommand("resubmit","JobType.maxMemoryMB=8000","--jobids", jobIds,dir=d)		
	else:
	    print("Status incomplete")
	    if res.get("statusFailureMsg"):
		print("Found reason for error: %s" % res.get("statusFailureMsg"))
