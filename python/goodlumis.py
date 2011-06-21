import copy
from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

# These numbers dictate how the rereco, prompt, DCS-only jsons are
# combined below.
last_rereco_run = 163869
last_prompt_run = 166861
assert last_prompt_run > last_rereco_run

# Lumis to manually throw out.
to_remove = {
    '166530': [[1,105]], # During physics-declared, this run had Mu30_v3 prescaled by 20. https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/PrescaleChanges?RUN=166530   https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/LumiSections?RUN=166530
    }
to_remove = LumiList(compactList=to_remove)

runs_to_remove_from_dcsonly = range(160404, last_prompt_run+1)
# These runs are <= last_prompt_run, but they were not actually
# considered in the certification for the latest prompt JSON. So,
# don't drop them from the DCS-only list when combining later.
#runs_to_remove_from_dcsonly.remove()

DCSOnly_ll           = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll.removeRuns(runs_to_remove_from_dcsonly)

Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-%i_7TeV_PromptReco_Collisions11_JSON.txt'          % last_prompt_run)
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-%i_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt' % last_prompt_run)

May10_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-%i_7TeV_May10ReReco_Collisions11_JSON.txt'          % last_rereco_run)
May10MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-%i_7TeV_May10ReReco_Collisions11_JSON_MuonPhys.txt' % last_rereco_run)

def combine(may10_ll, prompt_ll, dcsonly_ll=None):
    prompt_ll = copy.deepcopy(prompt_ll)
    prompt_ll.removeRuns(xrange(160404, last_rereco_run+1))
    ll = may10_ll | prompt_ll
    if dcsonly_ll is not None:
        dcsonly_ll = copy.deepcopy(dcsonly_ll)
        dcsonly_ll.removeRuns(runs_to_remove_from_dcsonly)
        ll = ll | dcsonly_ll
    return ll

Run2011A_ll          = combine(May10_ll,          Prompt_ll)
Run2011AMuonsOnly_ll = combine(May10MuonsOnly_ll, PromptMuonsOnly_ll)

Run2011APlusDCSOnly_ll          = combine(May10_ll,          Prompt_ll,          DCSOnly_ll)
Run2011APlusDCSOnlyMuonsOnly_ll = combine(May10MuonsOnly_ll, PromptMuonsOnly_ll, DCSOnly_ll)

Run2010_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON.txt')
Run2010MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_MuonPhys.txt')
                                   
all_ll_names = ['DCSOnly', 'DCSOnlyForNewRuns', 'Prompt', 'PromptMuonsOnly', 'May10', 'May10MuonsOnly', 'Run2011A', 'Run2011AMuonsOnly', 'Run2011APlusDCSOnly', 'Run2011APlusDCSOnlyMuonsOnly', 'Run2010', 'Run2010MuonsOnly']

def all_lls():
    return [(x, eval(x + '_ll')) for x in all_ll_names]

for base_name, ll in all_lls():
    exec '%s_ll = ll - to_remove' % base_name
    exec '%s = for_cmssw(%s_ll)' % (base_name, base_name)

if __name__ == '__main__':
    import sys
    if 'write' in sys.argv:
        Run2011APlusDCSOnlyMuonsOnly_ll.writeJSON('Run2011APlusDCSOnlyMuonsOnly.json')
    elif 'write_all' in sys.argv:
        for base_name, ll in all_lls():
            ll.writeJSON('zp2mu_goodlumis_%s.json' % base_name)
