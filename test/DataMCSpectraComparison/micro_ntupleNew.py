#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

path = 'data/Run2015MuonsOnly/ana_datamc_data.root'
#path = './zp2mu_histos.root'
tmp_fn = 'micro_ntupleNew.temp.txt'
#branch_spec = 'vertex_m'
branch_spec = 'run:lumi:event:vertex_m:'#dil_mass
#cut = 'OurSel2012'
#cut = 'lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 45 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 && (lep_triggerMatchPt[0] > 50 || lep_triggerMatchPt[1] > 50) && cos_angle > -0.9998   && vertex_chi2 < 10 && GoodData && OppSign && vertex_m>50'
cut = 'lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 53 && fabs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] > 0 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && (lep_pt_err[0] / lep_pt[0]) < 0.3 && (lep_sumPt[0] / lep_tk_pt[0]) < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 53 && fabs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] > 0 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && (lep_pt_err[1] / lep_pt[1]) < 0.3 && (lep_sumPt[1] / lep_tk_pt[1]) < 0.1 && ((lep_triggerMatchPt[0] > 50) || (lep_triggerMatchPt[1] > 50)) && cos_angle > -0.9998   && vertex_chi2 < 20 && GoodData && OppSign && vertex_m>50'
#cut = 'OurSel2012'

f = ROOT.TFile(path)
t = f.SimpleNtuplerstartup.Get('t')
t.GetPlayer().SetScanRedirect(True)
t.GetPlayer().SetScanFileName(tmp_fn)
#- colsize=10 is needed in order not to currupt the event number
t.Scan(branch_spec, cut, "colsize=10")
#t.Scan('*', cut)
t.GetPlayer().SetScanRedirect(False)
f.Close()

lines = [line.split(' *')[1:] for line in open(tmp_fn).readlines() if ' * ' in line and 'Row' not in line]
#- This loop is to keep only the first (highest-rank) dimuon in an event
prev_run = prev_lumi = prev_event = 0
cleaned_lines = []
for x in lines:
    cleaned_line = []
    curr_run   = x[0]
    curr_lumi  = x[1]
    curr_event = x[2]
    mass = x[3]
    #    print curr_run, curr_lumi, curr_event, mass
    if curr_run != prev_run or curr_lumi != prev_lumi or curr_event != prev_event:
        cleaned_line.append(mass)
        cleaned_lines.append(cleaned_line)
        prev_run   = curr_run
        prev_lumi  = curr_lumi
        prev_event = curr_event
    else:
        print 'delete line: ',curr_run,' ',curr_lumi,' ',curr_event,' mass ',mass
lines = ['\t'.join(y.strip() for y in x) for x in cleaned_lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
t = ROOT.TTree('t','')
t.ReadFile(tmp_fn, 'vertex_m')
f.Write()
f.Close()

#os.remove(tmp_fn)
