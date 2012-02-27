import gzip
from collections import defaultdict
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import to_pickle, from_pickle
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

def build_slots(lumi_per_slot): 
    slots = []
    lookup = {}
    
    curr_slot = defaultdict(list)
    tot_lumi_in_slot = 0

    curr_run = 0
    curr_ls = 0

    # lumiCalc2.py -i Run2011MuonsOnly.json -o Run2011MuonsOnly.lumibyls.csv lumibyls ; gzip Run2011MuonsOnly.lumibyls.csv
    for line in gzip.open('Run2011MuonsOnly.lumibyls.csv.gz'):
        try:
            line = line.split(',')     #Run,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub)    try:
            lumi = float(line[-1])/1e6
        except ValueError:
            print 'build_slots: skipping csv line:', line
            continue

        run = int(line[0])

        # the ls string will be of the format 'x:x' if the LS is selected
        # if not it might be 'x:0' or 'x:0:0' for some reason; recorded lumi will be 0 too
        ls = line[1].split(':')
        if len(ls) != 2 or ls[0] != ls[1]:
            assert int(lumi) == 0
            continue
        ls = int(ls[0])

        # ensure contiguity in slots. the csv from lumicalc should
        # already be in increasing (run,ls) order but check anyway.
        if run > curr_run:
            curr_run = run
            curr_ls = 0
        elif run < curr_run:
            raise ValueError('input not in increasing run order')
        if ls > curr_ls:
            curr_ls = ls
        elif ls < curr_ls:
            raise ValueError('input not in increasing ls order')
        
        curr_slot[run].append(ls)
        lookup[(run,ls)] = len(slots)

        tot_lumi_in_slot += lumi

        # if we've accumulated enough lumi, save this slot and make a new one
        if tot_lumi_in_slot > lumi_per_slot:
            slots.append((tot_lumi_in_slot, curr_slot))
            curr_slot = defaultdict(list)
            tot_lumi_in_slot = 0
        
    slots.append((tot_lumi_in_slot, curr_slot))
    print 'last slot has lumi', tot_lumi_in_slot
    return slots, lookup

lumi_per_slot = 491.7
slots, lookup = build_slots(lumi_per_slot)
to_pickle((slots, lookup), 'Run2011MuonsOnly.slots.%iipb.gzpickle' % lumi_per_slot)
#slots, lookup = from_pickle('Run2011MuonsOnly.slots.%iipb.gzpickle' % lumi_per_slot)
num_slots = len(slots)

mass_ranges = [
    ( 60,120),
    (120,200),
    (200,400),
    (120,600),
    ]

histos = [ROOT.TH1F('%s_%s' % m, 'event count for mass range %s-%s GeV' % m, num_slots, 0, num_slots) for m in mass_ranges]
    
f = ROOT.TFile('data/Run2011MuonsOnly/ana_datamc_data.root')
t = f.SimpleNtupler.Get('t')
#cut = ROOT.TTreeFormula('cut', 'OurSelNewNoSign && SameSign', t)
cut = ROOT.TTreeFormula('cut', 'OurSelNew', t)

for jentry, tt in ttree_iterator(t):
    if not cut.EvalInstance():
        continue

    slot = lookup[(t.run, t.lumi)]
    for i,(l,h) in enumerate(mass_ranges):
        if l <= t.dil_mass < h:
            histos[i].Fill(slot)

ROOT.gStyle.SetOptStat(10)
ps = plot_saver('plots/lumislots', log=False)
for h in histos:
    h.SetMarkerSize(1)
    h.SetMarkerStyle(0)
    h.SetMinimum(0)
    h.GetXaxis().SetTitle('lumi slot #')
    h.GetYaxis().SetTitle('events/%i pb^{-1}' % lumi_per_slot)
    h.Draw('hist e1 text90')
    ps.c.Update()
    s = h.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.73)
    s.SetY1NDC(0.92)
    s.SetX2NDC(0.98)
    s.SetY2NDC(0.99)
                    
    ps.save(h.GetName())
