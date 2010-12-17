import os
from FWCore.PythonUtilities.LumiList import LumiList

Run2010AB_ll          = Nov4_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON.txt')
Run2010ABMuonsOnly_ll = Nov4MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON_MuonPhys.txt')

Run2010AB = Run2010AB_ll.getCMSSWString().split(',')
Run2010ABMuonsOnly = Run2010ABMuonsOnly_ll.getCMSSWString().split(',')


