#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

path = 'data/Run2012MuonsOnly/ana_datamc_data.root'
tmp_fn = 'micro_ntuple.temp.txt'
branch_spec = 'vertex_m'
cut = 'OurSel2012'

f = ROOT.TFile(path)
t = f.SimpleNtupler.Get('t')
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

os.remove(tmp_fn)
