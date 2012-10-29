import copy
from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

# These run numbers guide the combination of the prompt and DCS-only
# JSONs.
first_run = 190456
last_rereco_run = 198523
last_prompt_run = 205618
last_run = 210000

# Sometimes the same run-range json gets made in other versions.
prompt_version = ''

# Lumis to manually throw out.
to_remove = {'190949': [[82,1149]], '191090': [[56,339]]}   # These are 20/pb of "low pileup" runs in which they enabled only Mu15 and disabled Mu40 (set prescale to 0).
to_remove.update({'191367': [[1,289]], '191391': [[1,14]]}) # Runs with < 0.2/pb data where triggers were switched off (prescale set to 0).  Just in DCS-only JSON.
to_remove.update({'193112': [[54,235]], '193116': [[1,693]]}) # ~2/pb + ~5/pb of "low pileup" runs in which they enabled only Mu15 and disabled Mu40 (set prescale to 0).
to_remove = LumiList(compactList=to_remove)

# These runs are <= last_prompt_run, but they were not actually
# considered in the certification for the latest prompt JSON. So,
# don't drop them from the DCS-only list when combining later.
holes = []

# Order of priorities: rereco, then prompt reco, and then DCS.
runs_to_remove_from_prompt  = range(first_run, last_rereco_run+1)
runs_to_remove_from_dcsonly = range(first_run, last_prompt_run+1)
for hole in holes:
    print 'goodlumis warning: re-adding "hole" run %i from DCS-only list' % hole
    runs_to_remove_from_dcsonly.remove(hole)

DCSOnly_ll           = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll.removeRuns(runs_to_remove_from_dcsonly)

# Remove runs outside the range [first_run, last_run] since DCS-only
# list includes HI runs, etc.
for ll in (DCSOnly_ll, DCSOnlyForNewRuns_ll):
    ll.removeRuns(xrange(1, first_run))
    ll.removeRuns(xrange(last_run+1, 300000)) # dummy number

# July 13th reprocessing of 2012A and 2012B
Jul13_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt')
Jul13MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_MuonPhys_v3.txt')

# August 6th reprocessing of 5 runs of 2012A
Aug06_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt')
Aug06MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON_MuonPhys.txt')

# August 24th reprocessing of 2012C v1
Aug24_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-%i_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'          % last_rereco_run)
Aug24MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-%i_8TeV_24Aug2012ReReco_Collisions12_JSON_MuonPhys.txt' % last_rereco_run)

# Prompt reconstruction, 2012C v2 and 2012D
Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_%i-%i_8TeV_PromptReco_Collisions12_JSON%s.txt'          % (first_run, last_prompt_run, prompt_version))
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_%i-%i_8TeV_PromptReco_Collisions12_JSON_MuonPhys%s.txt' % (first_run, last_prompt_run, prompt_version))

# A few runs with Level-1 calorimeter problems, good for muon physics
NoL1TMuonsOnly_ll_1 = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_195378-195379_8TeV_PromptReco_Collisions12_JSON_MuonPhys_NoL1T.txt')
NoL1TMuonsOnly_ll_2 = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_196353-196501_8TeV_PromptReco_Collisions12_JSON_MuonPhys_NoL1T.txt')
NoL1TMuonsOnly_ll_3 = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_200601-200601_8TeV_PromptReco_Collisions12_JSON_MuonPhys_NoL1T.txt')
NoL1TMuonsOnly_ll = [NoL1TMuonsOnly_ll_1, NoL1TMuonsOnly_ll_2, NoL1TMuonsOnly_ll_3]

def combine(prompt_ll, rereco1_ll, rereco2_ll, rereco3_ll, dcsonly_ll=None):
    prompt_ll = copy.deepcopy(prompt_ll)
    prompt_ll.removeRuns(runs_to_remove_from_prompt)
    ll = prompt_ll | rereco1_ll | rereco2_ll | rereco3_ll
    if dcsonly_ll is not None:
        dcsonly_ll = copy.deepcopy(dcsonly_ll)
        dcsonly_ll.removeRuns(runs_to_remove_from_dcsonly)
        ll = ll | dcsonly_ll
    return ll

# Combine all lists
Run2012_ll          = combine(Prompt_ll,          Jul13_ll,          Aug06_ll,          Aug24_ll)
Run2012MuonsOnly_ll = combine(PromptMuonsOnly_ll, Jul13MuonsOnly_ll, Aug06MuonsOnly_ll, Aug24MuonsOnly_ll)
for x in NoL1TMuonsOnly_ll:
    Run2012MuonsOnly_ll = Run2012MuonsOnly_ll | x

Run2012PlusDCSOnly_ll          = combine(Prompt_ll,          Jul13_ll,          Aug06_ll,          Aug24_ll,          DCSOnly_ll)
Run2012PlusDCSOnlyMuonsOnly_ll = combine(PromptMuonsOnly_ll, Jul13MuonsOnly_ll, Aug06MuonsOnly_ll, Aug24MuonsOnly_ll, DCSOnly_ll)
for x in NoL1TMuonsOnly_ll:
    Run2012PlusDCSOnlyMuonsOnly_ll = Run2012PlusDCSOnlyMuonsOnly_ll | x

all_ll_names = ['DCSOnly', 'DCSOnlyForNewRuns', 'Jul13', 'Jul13MuonsOnly', 'Prompt', 'PromptMuonsOnly', 'Run2012', 'Run2012MuonsOnly', 'Run2012PlusDCSOnly', 'Run2012PlusDCSOnlyMuonsOnly']

def all_lls():
    return [(x, eval(x + '_ll')) for x in all_ll_names]

for base_name, ll in all_lls():
    exec '%s_ll = ll - to_remove' % base_name
    exec '%s = for_cmssw(%s_ll)' % (base_name, base_name)

if __name__ == '__main__':
    import sys
    if 'write' in sys.argv:
        Run2012MuonsOnly_ll.writeJSON('Run2012MuonsOnly.json')
    elif 'write_all' in sys.argv:
        for base_name, ll in all_lls():
            ll.writeJSON('zp2mu_goodlumis_%s.json' % base_name)
