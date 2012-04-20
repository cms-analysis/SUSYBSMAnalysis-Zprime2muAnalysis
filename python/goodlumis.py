import copy
from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

# These run numbers guide the combination of the prompt and DCS-only
# JSONs.
first_run = 190456
last_prompt_run = 191276
last_run = 191849

# Sometimes the same run-range json gets made in other versions.
prompt_version = ''

# Lumis to manually throw out.
to_remove = {'190949': [[82,1149]], '191090': [[56,339]]} # These are 20/pb of "low pileup" runs in which they enabled only Mu15 and disabled Mu40 (set prescale to 0).
to_remove = LumiList(compactList=to_remove)

# These runs are <= last_prompt_run, but they were not actually
# considered in the certification for the latest prompt JSON. So,
# don't drop them from the DCS-only list when combining later.
holes = []

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

Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_%i-%i_8TeV_PromptReco_Collisions12_JSON%s.txt'          % (first_run, last_prompt_run, prompt_version))
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_%i-%i_8TeV_PromptReco_Collisions12_JSON_MuonPhys%s.txt' % (first_run, last_prompt_run, prompt_version))

def combine(prompt_ll, dcsonly_ll=None):
    ll = copy.deepcopy(prompt_ll)
    if dcsonly_ll is not None:
        dcsonly_ll = copy.deepcopy(dcsonly_ll)
        dcsonly_ll.removeRuns(runs_to_remove_from_dcsonly)
        ll = ll | dcsonly_ll
    return ll

Run2012_ll          = combine(Prompt_ll)
Run2012MuonsOnly_ll = combine(PromptMuonsOnly_ll)

Run2012PlusDCSOnly_ll          = combine(Prompt_ll,          DCSOnly_ll)
Run2012PlusDCSOnlyMuonsOnly_ll = combine(PromptMuonsOnly_ll, DCSOnly_ll)

all_ll_names = ['DCSOnly', 'DCSOnlyForNewRuns', 'Prompt', 'PromptMuonsOnly', 'Run2012', 'Run2012MuonsOnly', 'Run2012PlusDCSOnly', 'Run2012PlusDCSOnlyMuonsOnly']

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
