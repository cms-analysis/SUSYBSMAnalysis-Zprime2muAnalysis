import sys
from sys import argv
import os
sys.path.append('cfgs/')
import subprocess
import threading, Queue, time
verbose = False
users = {
	#"jan":["srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/jschulte/SingleMuon/"]
	#"jan":["srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/jschulte/"]
#	"jan":["gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/ADD/"]
	#"jan":["gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/CRAB_UserFiles/ana_datamc_dy_4Jets_2017"]
#	"jan":["gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"]
"fan":["root://cmseos.fnal.gov//store/user/zhangfa/DYMC2017resNOMUO/CRAB_UserFiles/dileptonAna_resolution_dy_4Jets_2017_nomuo"]
}


def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '%' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()

def removeFile(source, verbose=False):
    if (verbose):
        print 'gfal-rm ' + source

    subprocess.call(['gfal-rm ' + source], shell=True,stdout=open(os.devnull, 'wb'))
    return

def copyFile(source, destination, user, verbose=False):
    if (verbose):
        print 'gfal-copy -f  %s' + source + ' file:///' + destination
    #print 'srmcp srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv1?SFN=' + source + ' file:///' + destination
#lcg-cp -v -b -D srmv2 SURL  file://local_file
    subprocess.call(['gfal-copy -f %s'%source + ' file:///' + destination],shell=True,stdout=open(os.devnull, 'wb'))
    #time.sleep(1)
    return

class FetchingThread(threading.Thread):
    def __init__(self, verbose=False, replaceFilenameMap={}):
        self.verbose = verbose
        self.replaceFilenameMap = replaceFilenameMap
        threading.Thread.__init__(self)

    def run(self):
        global theThreadId
        self.threadId = theThreadId
        theThreadId += 1
        if (self.verbose): print "(thread " + str(self.threadId) + "): Thread started"

        nErrors = 0
        #rejectString = ""

        while (theListFiles.empty() == False):
            nextFile = theListFiles.get()
            if (self.verbose):
                print "(thread " + str(self.threadId) + "): Getting file " + os.path.basename(nextFile)

            # sort out log files

            #change file name if it is in self.replaceFilenameMap
            try:
                copyFile(nextFile, os.getcwd()+"/"+outDir+"/"+nextFile.split("/")[-1] , self.verbose)
                #if (os.path.exists(os.getcwd()+"/"+outDir+"/"+nextFile.split("/")[-1] ) == True):
                #    removeFile(nextFile, self.verbose)
                #else:
                #    print "ERROR: Copying failed: " + nextFile
                #    nErrors += 1
            except KeyboardInterrupt:
                print "Copying aborted: %s" % nextFile

        if (nErrors != 0):
            print "(thread " + str(self.threadId) + "): There were " + str(nErrors) + " copying errors!"
        if (self.verbose): print "(thread " + str(self.threadId) + "): Thread finished"

class StatusInfoThread(threading.Thread):
    def __init__(self, nTotalFiles, delay=10, verbose=False):
        self.nTotalFiles = nTotalFiles
        self.delay = delay
        self.verbose = verbose
        threading.Thread.__init__(self)

    def run(self):
        import sys
        self.threadId = "Status"
        if (self.verbose): print "(thread " + str(self.threadId) + "): Thread started"

        while (theListFiles.empty() == False):
            nFilesLeft = theListFiles.qsize()
            percentage = 100 * (self.nTotalFiles - nFilesLeft) / self.nTotalFiles
            sys.stdout.write("\r\033[1;34mStatus: %i%% (%i files left)\033[m" % (percentage, nFilesLeft))
            sys.stdout.flush()
            time.sleep(self.delay)

        if (self.verbose): print "(thread " + str(self.threadId) + "): Thread finished"

def __getCommandOutput2(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        if int(err) == 256:
            print 'Path in %s does not exist!' % (command)
            return ' '
        else:
            raise RuntimeError, '%s failed w/ exit code %d' % (command, err)
    return data

def getPathList(path):

		
	command = 'gfal-ls ' + path

    	output = __getCommandOutput2(command).splitlines()
	result = []

	for name in output:
		if '_' in name:
			result.append(path+'/'+name)

	return result

def getFileList(path,result):

		
	command = 'gfal-ls ' + path

    	output = __getCommandOutput2(command).splitlines()
	print output
	if len(output) > 0:
		if output[0] == path:
			result = output
			return result
		
		for name in output:
			if not 'root' in name and not 'failed' in name and not "log" in name:
				result = getFileList(path+'/'+name,result)
			elif not "failed" in name and not "log" in name:
				result.append(path+'/'+name)

#		if not "root" in output:
#			for newPath in output:
#				if not "failed" in newPath:
#					result = getFileList(path+"/"+newPath,result)			
#		else:
#			for index, file in enumerate(output):
#				output[index] = path+"/"+file
#			result = result + output
	return result

def mergeExpected(configName,tag,config):

	path = "results_%s%s/fromCRAB"%(configName,tag)
	
	
        for massRange in config.massesExp:
                mass = massRange[1]
		nJobs = massRange[3]
                while mass <= massRange[2]:
			fileList = []
			for i in range(0,nJobs):
                        	fileName = path + "/expectedLimit_%s%s_%s_%d.root"%(configName,tag,mass,i+1)
				if os.path.isfile(fileName):
					fileList.append(fileName)
				i+=1
			print "merging mass point %d GeV, %d/%d files present"%(mass,len(fileList),nJobs)	
			command = ["hadd","-f","%s/higgsCombine%s.MarkovChainMC.mH%d.123456.root"%(path,configName+tag,mass)]
			command += fileList
			subprocess.call(command,stdout=open(os.devnull, 'wb'))
                        mass += massRange[0]

def renameObs(configName,tag,config):

	path = "results_%s%s/fromCRAB"%(configName,tag)
	
	i = 0	
        for massRange in config.masses:
                mass = massRange[1]
                while mass < massRange[2]:
                        fileName = path + "/observedLimit_%s%s_%d.root"%(configName,tag,i+1)
			print "renaming result %d to mass point %d GeV"%(i+1,mass)	
			command = ["cp","%s"%fileName,"%s/higgsCombine%s.MarkovChainMC.mH%d.root"%(path,configName+tag,mass)]
			subprocess.call(command,stdout=open(os.devnull, 'wb'))
                        mass += massRange[0]
			i += 1



def main():
        import argparse
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    	parser.add_argument("--merge",dest="merge", action="store_true", default=False, help='merge expected limits')
        args = parser.parse_args()
	
	global outDir
	

	
	path = users['fan'][0]
	print "Searching for results in %s"%path

	paths = getPathList(path)
	for path in paths:
		#if not ("ADD" in path): continue
		#if not ("ZPrime" in path or "DYTo" in path or "ADD" in path): continue
		#path = paths[0]
		#outDir = path.split('/')[-1]
		outDir = ""
		#if not os.path.exists(outDir):
    		#	os.makedirs(outDir)
		files =  getFileList(path,[])
		print files
		nFiles = len(files)
		if nFiles > 0:
			print "Found %d files, begin downloading"%nFiles

	    		global theThreadId
	    		theThreadId = 1
	    		global theListFiles
	    		theListFiles = Queue.Queue(0)
    
   	 		# init file name replacing map
	    		replaceFilenameMap = {
	        		"download_1": 1
	        		}
        
    			# get directory content and fill file queue
		   	for file in files:
	    		#srmTools.getDir(theStoragePath, verbose):
            		# ignore directories
            			theListFiles.put(file)

	    		# start status thread
	       		StatusInfoThread(theListFiles.qsize(), 5, verbose).start()

	    		# start threads to copy all found files
			threads = []
	    		for iThread in range(0, 1):
				
        			threads.append(FetchingThread(verbose, replaceFilenameMap))
				threads[iThread].start()
        			# do not start all threads at the same time
        			time.sleep(5)

	    		for iThread in range(0, 1):
				
				threads[iThread].join()
        

		print "done here"


main()
