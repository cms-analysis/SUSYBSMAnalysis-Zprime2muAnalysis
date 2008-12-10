#!/usr/bin/env python

import os, sys

#print 'usage: %s [zp2mu_histos.root [resolution.ps]]' % sys.argv[0]    

# Run ROOT in batch mode.
sys.argv.append('-b')
from ROOT import *
sys.argv.remove('-b')

from PSDrawer import PSDrawer

if len(sys.argv) > 1:
    rootFile = sys.argv[1]
else:
    rootFile = 'zp2mu_histos.root'

if len(sys.argv) > 2:
    outputFile = sys.argv[2]
else:
    outputFile = 'resolution.ps'

if not os.path.isfile(rootFile):
    sys.stderr.write('input file %s does not exist!\n' % rootFile)
    sys.exit(1)

print 'Creating %s from %s...' % (outputFile, rootFile)

f = TFile(rootFile)
histos = f.Zprime2muResolution
psd = PSDrawer(outputFile)

#gStyle.SetHistMinimumZero(True)

def draw_page(page_type, histo_base_name, page_title, draw_log=False, fit_gaus=False, is_efficiency=False):
    if page_type == 'all':
        div = (2,5)
        levels = xrange(9)
    elif page_type == 'offline':
        div = (2,3)
        levels = xrange(4,9)
    elif page_type == 'no_gen':
        div = (2,4)
        levels = xrange(1,9)
    elif page_type == 'no_trig':
        div = (2,3)
        levels = [0] + range(4,9)
    else:
        raise ValueError, 'page_type %s not recognized' % page_type
    
    pad = psd.new_page(page_title, div)
    for i, level in enumerate(levels):
        subpad = pad.cd(i+1)
        if draw_log: subpad.SetLogy(1)

        h = histos.Get('%s%i' % (histo_base_name, level))

        if type(h) != type(TProfile):
            opt = 'hist'
        else:
            opt = ''
        
        if h is not None and is_efficiency:
            hgen = histos.Get('%s0' % histo_base_name)
            if hgen is not None:
                h = h.Clone()
                # Divide, computing binomial errors.
                h.Divide(h, hgen, 1, 1, 'B')
                h.SetMinimum(0.50)
                h.SetMaximum(1.02)
                opt += ' e'
                subpad.SetGrid(1)
            else:
                h = None
                
        if h is not None:
            h.Draw(opt)
            if fit_gaus:
                h.Fit('gaus', 'Q')

def draw(name, opt=''):
    h = histos.Get(name)
    if h:
        h.Draw(opt)
    return h


pad = psd.new_page('Mother particle id of leptons', (1,2))
for i in xrange(2):
    pad.cd(i+1).SetLogy(1)
    draw('LeptonOrigin%i' % i)

pad = psd.new_page('Dilepton acceptance', (2,2))
pad.cd(1)
hd = draw('GenMassAllEvents', 'hist')
pad.cd(2)
hn = draw('GenMassInAccept', 'hist')
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
            h.SetTitle(hgen.GetTitle().replace('Gen mass, GN,', 'Trigger eff. (' + '), L%i' % rec))
            h.Draw('hist e')

no_gen_pages = [
    ('EffVsEta', 'Reconstruction efficiency vs. #eta',  False, False, True),
    ('EffVsPhi', 'Reconstruction efficiency vs. #phi',  False, False, True),
    ('EffVsPt',  'Reconstruction efficiency vs. pT',    False, False, True),

    ('LeptonEtaDiff',    'Lepton #eta resolution'),
    ('LeptonPhiDiff',    'Lepton #phi resolution'),
    ('LeptonPtDiff',     'Lepton #DeltapT'),
    ('LeptonPtRes',      'Lepton pT resolution',   False, True),
    ('LeptonPRes',       'Lepton p resolution',    False, True),
    ('LeptonInvPtRes',   'Lepton 1/pT resolution', False, True),
    ('LeptonInvPRes',    'Lepton 1/p resolution',  False, True),
    ]

offline_pages = [
    ('LeptonInvPtResVPtGen',    'Lepton 1/pT resolution vs. gen pT'),
    ('LeptonInvPResVPGen',      'Lepton 1/p resolution vs. gen p'),
    
    ('LeptonInvPtPull',         'Lepton 1/pT pulls',              False, True),
    ('LeptonInvPPull',          'Lepton 1/p pulls',               False, True),
    ('LeptonInvPtResBarrel',    'Lepton 1/pT resolution, barrel', False, True),
    ('LeptonInvPResBarrel',     'Lepton 1/p resolution, barrel',  False, True),
    ('LeptonInvPtPullBarrel',   'Lepton 1/pT pulls, barrel',      False, True),
    ('LeptonInvPPullBarrel',    'Lepton 1/p pulls, barrel',       False, True),
    ('LeptonInvPtResEndcap',    'Lepton 1/pT resolution, endcap', False, True),
    ('LeptonInvPResEndcap',     'Lepton 1/p resolution, endcap',  False, True),
    ('LeptonInvPtPullEndcap',   'Lepton 1/pT pulls, endcap',      False, True),
    ('LeptonInvPPullEndcap',    'Lepton 1/p pulls, endcap',       False, True),
    
    ('DileptonMassRes',         'Dilepton mass resolution (dil - gen dil)',  False, True),
    ('DileptonResMassRes',      'Dilepton mass resolution (dil - gen res)',  False, True),
    ('ResonanceMassRes',        'Resonance mass resolution (res - gen res)', False, True),
    
    ('DileptonMassResVMass',    'Dilepton mass resolution (dil - gen dil) vs. mass'),
    ('DileptonResMassResVMass', 'Dilepton mass resolution (dil - gen res) vs. mass'),
    ('ResonanceMassResVMass',   'Resonance mass resolution (res - gen res) vs. mass'),

    ('ChargeDiff', 'Lepton charge assignment', True),
    ]

for page_info in no_gen_pages:
    draw_page('no_gen', *page_info)
for page_info in offline_pages:
    draw_page('offline', *page_info)

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
