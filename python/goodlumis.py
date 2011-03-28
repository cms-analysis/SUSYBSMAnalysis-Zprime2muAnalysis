from FWCore.PythonUtilities.LumiList import LumiList

DCSOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt_160405-161312')
DCSOnly = DCSOnly_ll.getCMSSWString().split(',')
