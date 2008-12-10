#!/usr/bin/env python

# usage: DrawZprime2muResolution.py [zp2mu_histos.root [resolution.ps]]

import os, sys
from PSDrawer import PSDrawer
from ROOT import TFile
 
if len(sys.argv) > 1: rootFile = sys.argv[1]
else:                 rootFile = 'zp2mu_histos.root'

if len(sys.argv) > 2: outputFile = sys.argv[2]
else:                 outputFile = 'resolution.ps'

if not os.path.isfile(rootFile):
    sys.stderr.write('input file %s does not exist!\n' % rootFile)
    sys.exit(1)

print 'Creating %s from %s...' % (outputFile, rootFile)

f = TFile(rootFile)
histos = f.Zprime2muResolution
psd = PSDrawer(outputFile)

################################################################################

pad = psd.new_page('Mother particle id of leptons', (1,2))
for i in xrange(2):
    pad.cd(i+1).SetLogy(1)
    psd.draw_if(histos, 'LeptonOrigin%i' % i)

pad = psd.new_page('Dilepton acceptance', (2,2))
pad.cd(1)
hd = psd.draw_if(histos, 'GenMassAllEvents', 'hist')
pad.cd(2)
hn = psd.draw_if(histos, 'GenMassInAccept', 'hist')
h = hn.Clone('GenMassAcceptance')
h.SetTitle('Acceptance vs. mass')
h.Divide(hd)
pad.cd(3)
h.Draw()

pad = psd.new_page('Trigger efficiency', (3,4))
gen_hists = [histos.Get('TrigEffVsDilMass0%i' % i) for i in range(3)]
for i, h in enumerate(gen_hists):
    if h is not None:
        pad.cd(i+1).SetGrid(1)
        h.Draw('hist')
for rec in xrange(1,4):
    for i, hgen in enumerate(gen_hists):
        pad.cd(3*rec + i + 1).SetGrid(1)
        hnum = histos.Get('TrigEffVsDilMass%i%i' % (rec, i))
        if hnum is not None:
            h = hnum.Clone()
            h.Divide(hnum, hgen, 1, 1, 'B')
            h.SetMinimum(0.89)
            h.SetMaximum(1.01)
            s = hgen.GetTitle().replace('Gen mass, GN,', '')
            h.SetTitle('Trigger eff. (%s), L%i' % (s, rec))
            h.Draw('hist e')

for var, title in [('Eta', '#eta'), ('Phi', '#phi'), ('Pt', 'pT')]:
    pad = psd.new_page('Reconstruction efficiency vs. %s' % title, (2,4))
    hgen = histos.Get('EffVs%s0' % var)
    if hgen is not None:
        for rec in xrange(1,9):
            pad.cd(rec).SetGrid(i)
            hnum = histos.Get('EffVs%s%i' % (var, rec))
            if hnum is not None:
                h = hnum.Clone()
                h.Divide(hnum, hgen, 1, 1, 'B')
                h.SetMinimum(0.50)
                h.SetMaximum(1.01)
                h.Draw('hist e')

def draw_efficiency(num_name, den_name, min=0.7, max=1.02):
    hnum = histos.Get(num_name)
    hden = histos.Get(den_name)
    if None not in (hnum, hden):
        h = hnum.Clone()
        h.Divide(hnum, hden, 1, 1, 'B')
        h.SetMinimum(min)
        h.SetMaximum(max)
        h.Draw('hist e')
    
# Total - acceptance, trigger and "off-line" reconstruction - efficiency.
pad = psd.new_page('Dilepton total rec. efficiency', (2,3))
for rec in xrange(4,9):
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%i1' % rec, 'DilRecEffVsMass00', min=0.5)

# "Off-line" efficiency to find at least two muons, relative to events
# in acceptance and accepted by the L1/HLT triggers.
pad = psd.new_page('Two-lepton offline rec. efficiency', (2,3))
for rec in xrange(4,9):
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%i0' % rec, 'DilRecEffVsMass01')

# "Off-line" efficiency to find opposite-sign dimuon, relative to
# events in acceptance and accepted by the L1/HLT triggers, and with
# two muons passing cuts.
pad = psd.new_page('Dilepton w/ cuts offline rec. efficiency', (2,3))
for rec in xrange(4,9):
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%i1' % rec, 'DilRecEffVsMass%i2' % rec)

# "Off-line" efficiency to find opposite-sign dimuon, relative to
# events in acceptance and accepted by the L1/HLT triggers.
pad = psd.new_page('Dilepton offline rec. efficiency', (2,3))
for rec in xrange(4,9):
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%i1' % rec, 'DilRecEffVsMass01')

psd.rec_level_page(histos, 'no_gen', 'LeptonEtaDiff',  'Lepton #eta resolution')
psd.rec_level_page(histos, 'no_gen', 'LeptonPhiDiff',  'Lepton #phi resolution')
psd.rec_level_page(histos, 'no_gen', 'LeptonPtDiff',   'Lepton #DeltapT')
psd.rec_level_page(histos, 'no_gen', 'LeptonPtRes',    'Lepton pT resolution',   fit_gaus=True)
psd.rec_level_page(histos, 'no_gen', 'LeptonPRes',     'Lepton p resolution',    fit_gaus=True)

psd.rec_level_page(histos, 'no_gen',  'LeptonInvPtRes',        'Lepton 1/pT resolution',            fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtResVPtGen',  'Lepton 1/pT resolution vs. gen pT')
psd.rec_level_page(histos, 'offline', 'LeptonInvPtResBarrel',  'Lepton 1/pT resolution, barrel',    fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtResEndcap',  'Lepton 1/pT resolution, endcap',    fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPull',       'Lepton 1/pT pulls',                 fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPullBarrel', 'Lepton 1/pT pulls, barrel',         fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPullEndcap', 'Lepton 1/pT pulls, endcap',         fit_gaus=True)

psd.rec_level_page(histos, 'no_gen',  'LeptonInvPRes',        'Lepton 1/p resolution',           fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPResVPGen',   'Lepton 1/p resolution vs. gen p')
psd.rec_level_page(histos, 'offline', 'LeptonInvPPull',       'Lepton 1/p pulls',                fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPResBarrel',  'Lepton 1/p resolution, barrel',   fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPResEndcap',  'Lepton 1/p resolution, endcap',   fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPPullBarrel', 'Lepton 1/p pulls, barrel',        fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPPullEndcap', 'Lepton 1/p pulls, endcap',        fit_gaus=True)

psd.rec_level_page(histos, 'offline', 'DileptonMassRes',    'Dilepton mass resolution (dil - gen dil)',  fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'DileptonResMassRes', 'Dilepton mass resolution (dil - gen res)',  fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'ResonanceMassRes',   'Resonance mass resolution (res - gen res)', fit_gaus=True)

psd.rec_level_page(histos, 'offline', 'DileptonMassResVMass',    'Dilepton mass resolution (dil - gen dil) vs. mass')
psd.rec_level_page(histos, 'offline', 'DileptonResMassResVMass', 'Dilepton mass resolution (dil - gen res) vs. mass')
psd.rec_level_page(histos, 'offline', 'ResonanceMassResVMass',   'Resonance mass resolution (res - gen res) vs. mass')

psd.rec_level_page(histos, 'offline', 'ChargeDiff', 'Lepton charge assignment', log_scale=True)

pad = psd.new_page('Charge misassignment vs. 1/pT', (2,3))
for rec in xrange(4,9):
    pad.cd(rec-3)
    hwrong = histos.Get('ChargeWrongVInvPt%i' % rec)
    hright = histos.Get('ChargeRightVInvPt%i' % rec)
    if None not in (hwrong, hright):
        h = hwrong.Clone()
        h.Divide(hright)
        h.SetTitle(h.GetTitle().replace('wrong q', '(wrong q)/(right q)'))
        h.Draw()

psd.close()
f.Close()
