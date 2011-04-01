from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

def combine(hww, prompt):
    hww.removeRuns(xrange(160404,161216))
    return prompt | hww

DCSOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')

HWW2011_ll = LumiList(compactList={"160431": [[19, 244]], "160432": [[2, 79]], "160454": [[1, 41]], "160455": [[1, 33]], "160456": [[1, 106]], "160808": [[7, 80]], "160815": [[1, 95]], "160819": [[1, 17], [19, 118]], "160827": [[1, 55]], "160835": [[1, 85], [87, 496], [569, 584]], "160853": [[62, 65]], "160871": [[68, 236]], "160872": [[1, 77]], "160873": [[1, 201]], "160874": [[1, 132]], "160875": [[1, 285]], "160876": [[1, 17]], "160877": [[1, 12], [79, 83]], "160888": [[45, 423]], "160890": [[1, 395]], "160894": [[1, 220]], "160898": [[25, 53]], "160907": [[15, 220]], "160911": [[1, 301]], "160913": [[1, 26]], "160914": [[1, 75]], "160915": [[1, 382]], "160916": [[1, 90]], "160935": [[33, 259]], "160936": [[1, 60]], "160937": [[1, 209]], "160938": [[1, 130]], "160939": [[1, 148]], "160940": [[1, 90]], "160942": [[1, 26]], "160943": [[1, 54]], "160954": [[78, 84]], "160955": [[1, 228]], "160956": [[1, 73]], "160957": [[1, 953]], "160994": [[44, 45]], "160998": [[1, 288]], "161008": [[1, 285]], "161016": [[1, 327]], "161020": [[1, 22]], "161103": [[1, 111]], "161106": [[1, 40]], "161107": [[1, 40]], "161113": [[1, 45]], "161116": [[1, 20]], "161117": [[1, 43]], "161119": [[1, 350]], "161165": [[1, 100]], "161176": [[1, 31]], "161217": [[1, 788]], "161222": [[1, 136]], "161223": [[1, 409]], "161224": [[1, 2]], "161233": [[33, 63]], "161310": [[39, 127]], "161311": [[1, 704]], "161312": [[1, 1027]]})

Prompt_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-161216_7TeV_PromptReco_Collisions11_JSON.txt')
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-161216_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt')

Run2011A_ll = Prompt_ll # 161079-161352 all BAD due to ES for now, so don't bother combining
Run2011AMuonsOnly_ll = combine(HWW2011_ll, PromptMuonsOnly_ll)

for x in ['DCSOnly', 'HWW2011', 'Prompt', 'PromptMuonsOnly', 'Run2011A', 'Run2011AMuonsOnly']:
    exec '%s = for_cmssw(%s_ll)' % (x,x)

