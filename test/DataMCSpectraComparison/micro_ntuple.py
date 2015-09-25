#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

#path = 'data/DCSOnly_Mu50/ana_datamc_data.root'
path = 'data/Run2015MuonsOnly_50ns/ana_datamc_data.root'
#path = './zp2mu_histos.root'
#path = './zp2mu_histos_370.root'
tmp_fn = 'micro_ntuple.temp.txt'
#branch_spec = 'run:lumi:event:vertex_m:lep_pt[1]:lep_pt[0]:lep_sumPt[1]:lep_tk_pt[1]:lep_sumPt[0]:lep_tk_pt[0]'#dil_mass
branch_spec = 'run:lumi:event:vertex_m'
#branch_spec = 'run:lumi:event:dil_mass:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_eta[1]:triggerMatched:trigger_match_0:trigger_match_1:lep_tk_numberOfValidPixelHits[0]:lep_tk_numberOfValidPixelHits[1]:lep_glb_numberOfValidPixelHits[0]:lep_glb_numberOfValidPixelHits[1]:lep_triggerMatchEta[0]:lep_triggerMatchEta[1]'
#branch_spec = 'run:lumi:event:dil_mass:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_eta[1]:triggerMatched:trigger_match_0:trigger_match_1'
#cut = 'vertex_m>400'
#cut = 'event ==  56306745 && lep_isGlobalMuon[0] && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && lep_pt[1] > 45 && abs(lep_dB[0]) < 0.2 &&  abs(lep_dB[1]) < 0.2  && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 &&  GoodData && OppSign && vertex_m>50'
########## 2 TeV
#cut = 'event ==  56306745 && lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 45 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && ((lep_triggerMatchPt[0] > 45) || (lep_triggerMatchPt[1] > 45)) && cos_angle > -0.9998   && vertex_chi2 < 10 && GoodData && OppSign && vertex_m>50 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 '
#cut = 'lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 45 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 && ((lep_triggerMatchPt[0] > 45 && abs(lep_triggerMatchEta[0]) < 2.1) || (lep_triggerMatchPt[1] > 45 && abs(lep_triggerMatchEta[1]) < 2.1)) && cos_angle > -0.9998   && vertex_chi2 < 20 && GoodData && OppSign && vertex_m>50'
### OR ov trigger filters ### offline cut only pt>45
#cut = 'run==251643 && lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 50 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 50 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 && ((lep_triggerMatchPt[0] > 50 ) || (lep_triggerMatchPt[1] > 50 )) && cos_angle > -0.9998   && vertex_chi2 < 10 && GoodData && OppSign && vertex_m>50'
#&& vertex_m>60 && vertex_m<120
#cut = 'loose_2012_0 && loose_2012_1'
cut = 'OurSel2012 && vertex_m>60'

f = ROOT.TFile(path)
t = f.SimpleNtuplertunepnew.Get('t')
#t = f.SimpleNtuplerstartup.Get('t')
t.GetPlayer().SetScanRedirect(True)
t.GetPlayer().SetScanFileName(tmp_fn)
t.Scan(branch_spec, cut)
t.GetPlayer().SetScanRedirect(False)
f.Close()

lines = [line.split(' *')[1:] for line in open(tmp_fn).readlines() if ' * ' in line and 'Row' not in line]
lines = ['\t'.join(y.strip() for y in x) for x in lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
t = ROOT.TTree('t','')
t.ReadFile(tmp_fn, branch_spec)
f.Write()
f.Close()

#os.remove(tmp_fn)
