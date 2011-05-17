from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

DCSOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')

DCSOnlyForNewRuns_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll.removeRuns(xrange(160404, 163869+1))

Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt')
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt')

Run2011A_ll = Prompt_ll
Run2011AMuonsOnly_ll = PromptMuonsOnly_ll

for x in ['DCSOnly', 'DCSOnlyForNewRuns', 'Prompt', 'PromptMuonsOnly', 'Run2011A', 'Run2011AMuonsOnly']:
    exec '%s = for_cmssw(%s_ll)' % (x,x)

if __name__ == '__main__':
    import sys
    if 'WriteGoodLumiJSONs' in sys.argv:
        for name in ['DCSOnly_ll', 'DCSOnlyForNewRuns_ll', 'Prompt_ll', 'PromptMuonsOnly_ll', 'Run2011A_ll', 'Run2011AMuonsOnly_ll']:
            obj = eval(name)
            obj.writeJSON('zp2mu_goodlumis_%s.json' % name.replace('_ll', ''))
