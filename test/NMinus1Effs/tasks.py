#!/usr/bin/env python

import sys, os, glob
import das_client
import json
from itertools import combinations
from FWCore.PythonUtilities.LumiList import LumiList
from SUSYBSMAnalysis.Zprime2muAnalysis.hadd import hadd
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn

just_testing = 'testing' in sys.argv
if just_testing:
    sys.argv.remove('testing')
    
try:
    cmd, extra = sys.argv[1].lower(), sys.argv[2:]
except IndexError:
    print 'usage: utils.py command [extra]'
    sys.exit(1)

def do(cmd):
    print cmd
    ret = []
    if not just_testing:
        cmds = cmd.split('\n') if '\n' in cmd else [cmd]
        for cmd in cmds:
            if cmd != '' and not cmd.startswith('#'):
                ret.append(os.system(cmd))
    if len(ret) == 1:
        ret = ret[0]
    return ret

#latest_dataset = '/SingleMu/Run2012A-PromptReco-v1/AOD'
#latest_dataset = '/SingleMu/Run2012B-PromptReco-v1/AOD'
#latest_dataset = '/SingleMu/Run2012C-PromptReco-v1/AOD'
#latest_dataset = '/SingleMu/Run2012C-PromptReco-v2/AOD'
#latest_dataset = '/SingleMu/Run2012D-PromptReco-v1/AOD'
#latest_dataset = '/SingleMu/Run2015A-PromptReco-v1/AOD'
#latest_dataset = '/SingleMu/Run2015A-PromptReco-v1/AOD'
#latest_dataset = '/ExpressPhysics/Run2015B-Express-v1/FEVT'
latest_dataset = '/SingleMuon/Run2015B-PromptReco-v1/AOD'
#lumi_masks = ['Run2012PlusDCSOnlyMuonsOnly', 'Run2012MuonsOnly'] #, 'DCSOnly', 'Run2012']
#lumi_masks = ['DCSOnly', 'Run2015', 'Run2015MuonsOnly'] #, 'Run2012PlusDCSOnlyMuonsOnly', 'DCSOnly']
#lumi_masks = ['DCSOnly','Run2015', 'Run2015MuonsOnly'] #, 'Run2012PlusDCSOnlyMuonsOnly', 'DCSOnly']
lumi_masks = ['Run2016MuonsOnly'] #, 'Run2012PlusDCSOnlyMuonsOnly', 'DCSOnly']


if cmd == 'setdirs':
    crab_dirs_location = extra[0]
    do('mkdir -p ' + os.path.join(crab_dirs_location, 'psets'))
    do('ln -s %s crab' % crab_dirs_location)
    do('ln -s . crab/crab') # this is so crab -publish won't complain about being unable to find the pset if you launch it from crab/
    do('mkdir crab/publish_logs')

elif cmd == 'maketagdirs':
    extra = extra[0]
    do('rm data mc plots')
    for which in ['data', 'mc', 'plots']:
        #d = '~/nobackup/zp2mu_ana_datamc_%s/%s' % (which,extra)
        d = './zp2mu_ana_nminus1_%s/%s' % (which,extra)
        do('mkdir -p %s' % d)
        do('ln -s %s %s' % (d, which))

elif cmd == 'checkevents':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('grep TrigReport crab/crab_nminus1_%s/res/*stdout | grep \' p$\' | sed -e "s/ +/ /g" | awk \'{ s += $4; t += $5; u += $6; } END { print "summary: total: ", s, "passed: ", t, "failed: ", u }\'' % sample.name)

elif cmd == 'checkstatus':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab status -d crab/crab_ana_nminus1_%(name)s ' % sample)
        
elif cmd == 'report':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab report -d crab/crab_ana_nminus1_%(name)s ' % sample)
        
elif cmd == 'resubmit':
     from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
     for sample in samples:
         print sample.name
         do('crab resubmit -d crab/crab_ana_nminus1_%(name)s ' % sample)

elif cmd == 'kill':
     from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
     for sample in samples:
         print sample.name
         do('crab kill -d crab/crab_ana_nminus1_%(name)s ' % sample)

elif cmd == 'getoutput':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        print sample.name
        do('crab getoutput -d crab/crab_ana_nminus1_%(name)s --checksum=no' % sample)

#elif cmd == 'publishmc':
#    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
#    for sample in samples:
#        do('crab -c crab/crab_nminus1_%(name)s -publish >& crab/publish_logs/publish.crab_nminus1_%(name)s &' % sample)

elif cmd == 'anadatasets':
    print 'paste this into python/MCSamples.py:\n'
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    for sample in samples:
        ana_dataset = None
        fn = 'crab/publish_logs/publish.crab_nminus1_%s' % sample.name
        # yay fragile parsing
        for line in open(fn):
            if line.startswith(' total events'):
                ana_dataset = line.split(' ')[-1].strip()
                break
        if ana_dataset is None:
            raise ValueError('could not find ana_dataset from %s' % fn)
        print '%s.ana_dataset = "%s"' % (sample.name, ana_dataset)

elif cmd == 'gathermc':
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples
    extra = '_' + extra[0] if extra else ''
    for sample in samples:
        name = sample.name
        pattern = 'crab/crab_ana%(extra)s_nminus1_%(name)s/results/zp2mu_histos*root' % locals()
        fn = 'ana_nminus1_%(name)s.root' % locals()
        n = len(glob.glob(pattern))
        if n == 0:
            big_warn('no files matching %s' % pattern)
        else:
            files = glob.glob('crab/crab_ana%(extra)s_nminus1_%(name)s/results/zp2mu_histos*root' % locals())
            hadd('mc/ana_nminus1_%s.root' % name, files)

elif cmd == 'gatherdata':
    extra = (extra[0] + '_') if extra else ''

    for lumi_mask in lumi_masks:
        print lumi_mask
#        dirs = glob.glob('crab/crab_ana_datamc_%s_ExpressPhysicsRun2015B*' % lumi_mask)
        dirs = glob.glob('crab/crab_ana_nminus1_SingleMuonRun2016*' )
#        dirs = glob.glob('crab/crab_ana_datamc_%s_ExpressPhysicsRun2015B-Express_251161_251252' % lumi_mask)
        files = []
        for d in dirs:
            files += glob.glob(os.path.join(d, 'results/*.root'))

        wdir = os.path.join('data', lumi_mask)
        os.mkdir(wdir)
        hadd(os.path.join(wdir, 'ana_nminus1_data.root'), files)

#         for dir in dirs:
#             do('crab status -d %(dir)s ; crab report -d %(dir)s ' % locals())

        jsons = [os.path.join(dir, 'results/processedLumis.json') for dir in dirs]
        print jsons
        lls = [(j, LumiList(j)) for j in jsons]
        for (j1, ll1), (j2, ll2) in combinations(lls, 2):
            cl = (ll1 & ll2).getCompactList()
            print 'checking overlap between', j1, j2,
            if cl:
                raise RuntimeError('\noverlap between %s and %s lumisections' % (j1,j2))
            else:
                print cl
        
        reduce(lambda x,y: x|y, (LumiList(j) for j in jsons)).writeJSON('%(wdir)s/ana_nminus1_data.forlumi.json' % locals())
        do('brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /pb  -i %(wdir)s/ana_nminus1_data.forlumi.json  > %(wdir)s/ana_nminus1_data.lumi -b "STABLE BEAMS" ' % locals())
        do('tail -5 %(wdir)s/ana_nminus1_data.lumi' % locals())
        print 'done with', lumi_mask, '\n'
		#do('lumiCalc2.py -i %(wdir)s/ana_nminus1_data.forlumi.json overview > %(wdir)s/ana_nminus1_data.lumi' % locals())
		#do('pixelLumiCalc.py -i %(wdir)s/ana_nminus1_data.forlumi.json overview > %(wdir)s/ana_nminus1_data.lumi' % locals())
		#do('brilcalc lumi -i %(wdir)s/ana_nminus1_data.forlumi.json -n 0.962 > %(wdir)s/ana_nminus1_data.lumi' % locals())
		

		
		   

elif cmd == 'runrange':
    #cmd = 'dbs search --query="find min(run),max(run) where dataset=%s"' % latest_dataset
    #do(cmd)
    cmd = 'das_client.py --limit=0 --query="run dataset=%s"' % latest_dataset
    #do(cmd)
    runs=[]
    for line in os.popen(cmd):
        runs.append(int(line))
    #print runs
    print min(runs), ' ',max(runs)

    #query = "run dataset=%s" % latest_dataset
    #print query
    #data = das_client.json.loads(das_client.get_data('https://cmsweb.cern.ch',query,0,0,0))

elif cmd == 'checkavail':
    #cmd = 'dbs search --query="find run,lumi where dataset=%s"' % latest_dataset
    cmd = 'das_client.py --limit=0 --query="run,lumi dataset=%s"' % latest_dataset
    print '\n',cmd

    lumis = {}
    query = "run,lumi dataset="+latest_dataset
    jsondict=das_client.get_data('https://cmsweb.cern.ch',query,0,0,0)
    status = jsondict['status']
    if status != 'ok':
        print "DAS query status: %s"%(status)
    data = jsondict['data']
    for aRaw in data:
        run=aRaw['run'][0]['run_number']
        run_lumis=aRaw['lumi'][0]['number']
        if(len(run_lumis)):
#        for r in aRaw['run']:
#            run_2=r['run_number']
#        for l in aRaw['lumi']:
#            lumis_2=l['number']
            #print run, run_lumis
            #lumis.append((run,run_lumis))
            lumis[run]=run_lumis
#    print lumis

    ll = LumiList(compactList=lumis)
#print "ll", ll
    runrange = sorted(int(x) for x in ll.getCompactList().keys())

    #dcs_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY.txt') # JMTBAD import from goodlumis
    dcs_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt') # JMTBAD import from goodlumis

    #print "dcs_ll", dcs_ll
    dcs_runrange = sorted(int(x) for x in dcs_ll.getCompactList().keys())

    dcs_ll.removeRuns(xrange(dcs_runrange[0], runrange[0]))
    dcs_ll.removeRuns(xrange(runrange[-1]+1, dcs_runrange[-1]))

#- 2012A:    190450-193686
#- 2012B:    193752-196531
#- 2012C v1: 197556-198913
#- 2012C v2: 198934-203772
#- 2012D:    203773-
#- Reasons for excluding some runs:
#    191350, 192989, 192890, 193091, 204900, 206251, 207871: VdM scan
#    193092: very low pile-up run
#ok = LumiList(compactList={ "208540": [[99, 101]], "208551": [[581, 586]]})#fake list
    ok = LumiList(compactList={
#              "246958": [[1, 61]],
#              "246959": [[1, 77]],
#              "246960": [[1, 1], [37, 44]],
#              "246963": [[1, 11], [13, 13], [15, 312], [314, 314], [316, 448]],
#              "246965": [[1, 5]],
#              "247047": [[1, 5], [13, 13], [15, 15]],
#              "247049": [[3, 5]],
#              "247052": [[1, 60]],
#              "247054": [[2, 5], [14, 79]],
#              "247056": [[1, 5], [13, 18]],
#              "247057": [[1, 5], [13, 18], [26, 26]],
#              "247059": [[1, 3]],
#              "247060": [[1, 2], [4, 4]],
#              "247063": [[1, 5], [13, 18]],
#              "247068": [[1, 104], [106, 106], [108, 116], [118, 118], [121, 133]]
                  })#from lumiSummary.json

    print 'run range for', latest_dataset, ':', runrange[0], runrange[-1]
    print 'these lumis are in the DCS-only JSON but not (yet) in', latest_dataset
    print str(dcs_ll - ll - ok)

#elif cmd == 'oklist':  not true
#    addOk_ll = LumiList('crab/crab_datamc_EGammaRun2015A-Prompt_246958_247068_20150622005921/results/lumiSummary.json')
#    print addOk_ll

elif cmd == 'drawall':
    extra = extra[0] if extra else ''
    for lumi_mask in lumi_masks:
        r = do('python draw.py data/ana_nminus1_%s %s > out.draw.%s' % (lumi_mask,extra,lumi_mask))
        if r != 0:
            sys.exit(r)
    do('mv out.draw.* plots/')
    do('tlock ~/asdf/plots.tgz plots/nminus1_* plots/out.draw.*')

else:
    raise ValueError('command %s not recognized!' % cmd)
