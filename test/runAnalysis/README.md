README for the new user interface to run jobs in this framework

Currently, running jobs for the HistosFromPAT histogrammer (+ SimpleNtupler if desired) and the resolution studies (ResolutonUsingMC and ResolutionAtZ modules) is supported. 

Everything is steered through runAnalyis.py

The tool supports three basic modes of operation:
	-w option will create a CMSSW config (cmssw_cfg.py), but will not execute anything
	-l option will execute said file locally, running over the file(s) specified in setup.py 
	-s option will submit jobs to crab in a loop over a giving set of samples

By default, the histogrammer will be run, and nTuples can be added using the --addNTuples option. To change to a different tool, the available options are:
	-r produces histograms for resolution studies
	- ....

Further options:
	- By default, 2017 data or MC are processed. This can be changed using the --2016 or --2018 options
	- By default, the tool is run configured for dimuon pairs. It can be switched to use electrons using the -e option
	- By default, MC is process, use the -d option to switch to data
	- If MC is processed, background samples are used by default. For signal samples. To process signal, additional options have be set. --ci2016 and --ci2017 are supported right now for the non-resonant analysis
	
The available MC samples are collected as lists in samples.py 

Note that the appropriate generator filter are added to the cmssw config based on the sample name when running runAnalysis.py 

To check on crab jobs, resubmit failed jobs and harvest and merge output of succesful jobs, just run crabController.py
 
