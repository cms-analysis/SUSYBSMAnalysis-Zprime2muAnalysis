import copy
from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

DCSOnly_ll           = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll.removeRuns(xrange(160404, 165121+1))

Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-165121_7TeV_PromptReco_Collisions11_JSON.txt')
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-165121_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt')

May10_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON.txt')
May10MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys.txt')

def combine(may10_ll, prompt_ll):
    prompt_ll = copy.deepcopy(prompt_ll)
    prompt_ll.removeRuns(xrange(160404, 163869+1))
    return may10_ll | prompt_ll

Run2011A_ll          = combine(May10_ll,          Prompt_ll)
Run2011AMuonsOnly_ll = combine(May10MuonsOnly_ll, PromptMuonsOnly_ll)

Run2010_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON.txt')
Run2010MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_MuonPhys.txt')

for x in ['DCSOnly', 'DCSOnlyForNewRuns', 'Prompt', 'PromptMuonsOnly', 'May10', 'May10MuonsOnly', 'Run2011A', 'Run2011AMuonsOnly', 'Run2010', 'Run2010MuonsOnly']:
    exec '%s = for_cmssw(%s_ll)' % (x,x)

if __name__ == '__main__':
    import sys
    if 'WriteGoodLumiJSONs' in sys.argv:
        for name in ['DCSOnly_ll', 'DCSOnlyForNewRuns_ll', 'Prompt_ll', 'PromptMuonsOnly_ll', 'May10', 'May10MuonsOnly', 'Run2011A_ll', 'Run2011AMuonsOnly_ll']:
            obj = eval(name)
            obj.writeJSON('zp2mu_goodlumis_%s.json' % name.replace('_ll', ''))
