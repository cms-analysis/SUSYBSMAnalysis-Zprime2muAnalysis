#!/usr/bin/env python

import sys, os, glob
from FWCore.PythonUtilities.LumiList import LumiList
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
    if not just_testing:
        cmds = cmd.split('\n') if '\n' in cmd else [cmd]
        for c in cmds:
            if cmd != '' and not cmd.startswith('#'):
                os.system(c)

latest_dataset = '/SingleMu/Run2011A-PromptReco-v4/AOD'

if cmd == 'setdirs':
    do('''ln -s /uscms/home/tucker/nobackup/crab_dirs crab
mkdir -p psets
ln -s `pwd`/psets crab/psets''')

elif cmd == 'checkevents':
    from samples import samples
    for sample in samples:
        print sample.name
        do('grep TrigReport crab/crab_datamc_%s/res/*stdout | grep \' p$\' | sed -e "s/ +/ /g" | awk \'{ s += $4; t += $5; u += $6; } END { print "summary: total: ", s, "passed: ", t, "failed: ", u }\'' % sample.name)

elif cmd == 'publishmc':
    from samples import samples
    for sample in samples:
        do('crab -c crab/crab_datamc_${x} -publish >&! crab/publish_logs/publish.${x} &')

elif cmd == 'hadd':
    from samples import samples
    extra = '_' + extra[0] if extra else ''
    for sample in samples:
        name = sample.name
        pattern = 'crab/crab_ana%(extra)s_datamc_%(name)s/res/zp2mu_histos*root' % locals()
        fn = 'ana_datamc_%(name)s.root' % locals()
        n = len(glob.glob(pattern))
        if n == 0:
            big_warn('no files matching %s' % pattern)
        elif n == 1:
            do('cp %s %s' % (pattern, fn))
        else:
            do('hadd ana_datamc_%(name)s.root crab/crab_ana%(extra)s_datamc_%(name)s/res/zp2mu_histos*root' % locals())

elif cmd == 'gatherhistos':
    extra = extra[0] if extra else 'renameme'
    which = 'Run2011APlusDCSOnlyMuonsOnly'
    #which = 'Run2011A'
    dirs = glob.glob('crab/crab_ana_datamc_%s_SingleMu2011A_*' % which)
    files_glob = ' '.join([os.path.join(x, 'res/*.root') for x in dirs])
    do('''
mkdir -p ana_datamc_%(extra)s
ln -s /uscms_data/d2/tucker/zp2mu_ana_datamc_mc/latest ana_datamc_%(extra)s/mc
hadd ana_datamc_%(extra)s/ana_datamc_data.root %(files_glob)s
''' % locals())
    for dir in dirs:
        do('crab -c %s -report' % dir)
    jsons = [os.path.join(dir, 'res/lumiSummary.json') for dir in dirs]
    for i,json1 in enumerate(jsons):
        for json2 in jsons[i+1:]:
            # probably a better way to do this using LumiList but oh well
            print 'checking overlap between', json1, json2
            if os.popen('compareJSON.py --and %s %s' % (json1, json2)).read() != '{}\n':
                raise RuntimeError('overlap between %s and %s lumisections' % (json1, json2))
    jsons = ' '.join(jsons)
    do('''
mergeJSON.py %(jsons)s --output ana_datamc_%(extra)s/ana_datamc_data.forlumi.json
lumiCalc.py -i ana_datamc_%(extra)s/ana_datamc_data.forlumi.json overview > ana_datamc_%(extra)s/ana_datamc_data.lumi
''' % locals())

elif cmd == 'mclinks':
    extra = extra[0] if extra else 'renameme'
    for x in ['muonsonly', 'allgood']:
        print('ln -s ~/nobackup/ana_datamc_mc/%s %s' % (extra, x))

elif cmd == 'runrange':
    cmd = 'dbs search --query="find min(run),max(run) where dataset=%s"' % latest_dataset
    do(cmd)

elif cmd == 'checkavail':
    cmd = 'dbs search --query="find run,lumi where dataset=%s"' % latest_dataset
    print '\n',cmd

    lumis = []
    for line in os.popen(cmd):
    #for line in open('temp'):
        line = line.strip().split()
        if len(line) != 2:
            continue
        try:
            a = int(line[0])
            b = int(line[1])
        except ValueError:
            continue
        lumis.append((a,b))

    ll = LumiList(lumis=lumis)
    runrange = sorted(int(x) for x in ll.getCompactList().keys())

    dcs_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
    dcs_runrange = sorted(int(x) for x in dcs_ll.getCompactList().keys())

    dcs_ll.removeRuns(xrange(dcs_runrange[0], runrange[0]))
    dcs_ll.removeRuns(xrange(runrange[-1]+1, dcs_runrange[-1]))

    ok = LumiList(compactList={"165098": [[188, 189], [194, 194], [249, 249], [255, 255], [332, 332], [368, 368], [416, 440]],
                               "165099": [[107, 126]],
                               "165120": [[99, 109]],
                               "165205": [[249, 260]],
                               "165364": [[112, 113], [808, 808], [1331, 1352]],
                               "165400": [[80, 88]],
                               "165415": [[641, 642], [708, 711], [778, 779]],
                               "165467": [[710, 736]],
                               "165472": [[185, 185]],
                               "165486": [[103, 118]],
                               "165487": [[152, 154]],
                               "165506": [[171, 171]],
                               "165514": [[568, 570]],
                               "165523": [[60, 92]],
                               "165525": [[33, 68]],
                               "165537": [[297, 304]],
                               "165542": [[176, 189]],
                               "165548": [[364, 364], [590, 630]],
                               "165567": [[110, 113], [310, 314], [632, 637]],
                               "165570": [[3, 4], [84, 87], [947, 947]],
                               "165617": [[53, 53], [144, 144], [289, 291]],
                               "165619": [[20, 21], [23, 23]],
                               "165620": [[25, 25], [59, 63]],
                               "165633": [[318, 318]],
                               "166010": [[174, 175], [179, 180]],
                               "166033": [[1234, 1244]],
                               "166034": [[363, 363]],
                               "166049": [[886, 887], [918, 918]],
                               "166149": [[3, 24]],
                               "166161": [[121, 121], [124, 125], [127, 128]],
                               "166163": [[13, 13], [35, 35]],
                               "166374": [[65, 65], [189, 193]],
                               "166377": [[57, 57]],
                               "166380": [[712, 714]],
                               "166486": [[96, 96]],
                               })

    print 'run range for', latest_dataset, ':', runrange[0], runrange[-1]
    print 'these lumis are in the DCS-only JSON but not (yet) in', latest_dataset
    print str(dcs_ll - ll - ok)
