#!/usr/bin/env python

# usage: DrawZprime2muResolution.py [zp2mu_histos.root [resolution.ps]]

import os, sys
from PSDrawer import PSDrawer
from ROOT import TFile, TLatex, gStyle
 
if len(sys.argv) > 1: rootFile = sys.argv[1]
else:                 rootFile = 'zp2mu_histos.root'

if len(sys.argv) > 2:
    outputFile = sys.argv[2]
else:
    outputFile = os.path.basename(rootFile).replace('.root', '.resolution.ps')

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

def draw_efficiency(num_name, den, min=0.7, max=1.02):
    hnum = histos.Get(num_name)
    if type(den) == type(''):
        hden = histos.Get(den)
    else:
        hden = den
        
    if None not in (hnum, hden):
        h = hnum.Clone()
        h.Divide(hnum, hden, 1, 1, 'B')
        h.SetMinimum(min)
        h.SetMaximum(max)
        h.Draw('hist e')

        # Calculate and draw total efficiency.
        nden = hden.GetEntries()
        if nden > 0:
            effpct = 100.0*hnum.GetEntries()/nden
            eff = '#varepsilon = %.1f%%' % effpct
            t = TLatex()
            t.SetTextSize(0.1)
            t.SetNDC()
            t.DrawLatex(0.5, 0.3, eff)

        return h

gStyle.SetOptStat(1111)

pad = psd.new_page('Trigger efficiency', (3,4))
gen_hists = [histos.Get('TrigEffVsDilMass0%i' % i) for i in range(3)]
for i, h in enumerate(gen_hists):
    if h is not None:
        pad.cd(i+1).SetGrid(1)
        h.Draw('hist')
for rec in xrange(psd.TRIG_START, psd.OFFLINE_START):
    for i, hgen in enumerate(gen_hists):
        pad.cd(3*rec + i + 1).SetGrid(1)
        h = draw_efficiency('TrigEffVsDilMass%i%i' % (rec, i), hgen, 0.89, 1.01)
        if h is not None:
            h.SetTitle('Trigger eff. (%s), L%i' % (hgen.GetTitle().replace('Gen mass, GN,', ''), rec))

div, levels = psd.div_levels('no_gen')
for var, title in [('Eta', '#eta'), ('Phi', '#phi'), ('Pt', 'pT')]:
    pad = psd.new_page('Reconstruction efficiency vs. %s' % title, div)
    hgen = histos.Get('EffVs%s0' % var)
    if hgen is not None:
        for rec in levels:
            pad.cd(rec).SetGrid(i)
            draw_efficiency('EffVs%s%X' % (var, rec), hgen, 0.5, 1.01)

div, levels = psd.div_levels('offline')

# Total - acceptance, trigger and "off-line" reconstruction - efficiency.
pad = psd.new_page('Dilepton total rec. efficiency', div)
for rec in levels:
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%X1' % rec, 'DilRecEffVsMass00', min=0.5)

# "Off-line" efficiency to find at least two muons, relative to events
# in acceptance and accepted by the L1/HLT triggers.
pad = psd.new_page('Two-lepton offline rec. efficiency', div)
for rec in levels:
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%X0' % rec, 'DilRecEffVsMass01')

# "Off-line" efficiency to find opposite-sign dimuon, relative to
# events in acceptance and accepted by the L1/HLT triggers, and with
# two muons passing cuts.
pad = psd.new_page('Dilepton w/ cuts offline rec. efficiency', div)
for rec in levels:
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%X1' % rec, 'DilRecEffVsMass%X2' % rec)

# "Off-line" efficiency to find opposite-sign dimuon, relative to
# events in acceptance and accepted by the L1/HLT triggers.
pad = psd.new_page('Dilepton offline rec. efficiency', div)
for rec in levels:
    pad.cd(rec-3).SetGrid(1)
    draw_efficiency('DilRecEffVsMass%X1' % rec, 'DilRecEffVsMass01')

gStyle.SetOptStat(111111)

psd.rec_level_page(histos, 'no_gen', 'LeptonEtaDiff',  'Lepton #eta resolution', fit_gaus=True)
psd.rec_level_page(histos, 'no_gen', 'LeptonPhiDiff',  'Lepton #phi resolution', fit_gaus=True)
psd.rec_level_page(histos, 'no_gen', 'LeptonPtDiff',   'Lepton #DeltapT')
psd.rec_level_page(histos, 'no_gen', 'LeptonPtRes',    'Lepton pT resolution',   fit_gaus=True)
psd.rec_level_page(histos, 'no_gen', 'LeptonPRes',     'Lepton p resolution',    fit_gaus=True)

psd.rec_level_page(histos, 'no_gen',  'LeptonInvPtRes',          'Lepton 1/pT resolution',            fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'RMSLeptonInvPtResVPtGen', 'Lepton 1/pT resolution vs. gen pT', hist_cmds=[('SetMaximum', 0.5)])
psd.rec_level_page(histos, 'offline', 'LeptonInvPtResBarrel',    'Lepton 1/pT resolution, barrel',    fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtResEndcap',    'Lepton 1/pT resolution, endcap',    fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPull',         'Lepton 1/pT pulls',                 fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPullBarrel',   'Lepton 1/pT pulls, barrel',         fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPtPullEndcap',   'Lepton 1/pT pulls, endcap',         fit_gaus=True)

psd.rec_level_page(histos, 'no_gen',  'LeptonInvPRes',         'Lepton 1/p resolution',           fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'RMSLeptonInvPResVPGen', 'Lepton 1/p resolution vs. gen p', hist_cmds=[('SetMaximum', 0.5)])
psd.rec_level_page(histos, 'offline', 'LeptonInvPPull',        'Lepton 1/p pulls',                fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPResBarrel',   'Lepton 1/p resolution, barrel',   fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPResEndcap',   'Lepton 1/p resolution, endcap',   fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPPullBarrel',  'Lepton 1/p pulls, barrel',        fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'LeptonInvPPullEndcap',  'Lepton 1/p pulls, endcap',        fit_gaus=True)

psd.rec_level_page(histos, 'offline', 'DileptonMassRes',    'Dilepton mass resolution (dil - gen dil)',  fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'DileptonResMassRes', 'Dilepton mass resolution (dil - gen res)',  fit_gaus=True)
psd.rec_level_page(histos, 'offline', 'ResonanceMassRes',   'Resonance mass resolution (res - gen res)', fit_gaus=True)

psd.rec_level_page(histos, 'offline', 'RMSDileptonMassResVMass',    'Dilepton mass resolution (dil - gen dil) vs. mass')
psd.rec_level_page(histos, 'offline', 'RMSDileptonResMassResVMass', 'Dilepton mass resolution (dil - gen res) vs. mass')
psd.rec_level_page(histos, 'offline', 'RMSResonanceMassResVMass',   'Resonance mass resolution (res - gen res) vs. mass')

psd.rec_level_page(histos, 'offline', 'ChargeDiff', 'Lepton charge assignment', log_scale=True)

pad = psd.new_page('Charge misassignment vs. 1/pT', div)
for rec in levels:
    pad.cd(rec-3)
    hwrong = histos.Get('ChargeWrongVInvPt%X' % rec)
    hright = histos.Get('ChargeRightVInvPt%X' % rec)
    if None not in (hwrong, hright):
        h = hwrong.Clone()
        h.Divide(hright)
        h.SetTitle(h.GetTitle().replace('wrong q', '(wrong q)/(right q)'))
        h.Draw()

psd.close()
f.Close()
