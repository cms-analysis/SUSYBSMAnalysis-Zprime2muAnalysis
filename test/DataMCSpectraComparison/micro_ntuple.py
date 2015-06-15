#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

path = 'data/Run2012MuonsOnly/ana_datamc_data.root'
tmp_fn = 'micro_ntuple.temp.txt'
branch_spec = 'run:lumi:event:vertex_m'
cut = 'OurSel2012'

f = ROOT.TFile(path)
t = f.SimpleNtupler.Get('t')
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
#    else:
#        print 'delete line'
lines = ['\t'.join(y.strip() for y in x) for x in cleaned_lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
t = ROOT.TTree('t','')
t.ReadFile(tmp_fn, 'vertex_m')
f.Write()
f.Close()

#os.remove(tmp_fn)
