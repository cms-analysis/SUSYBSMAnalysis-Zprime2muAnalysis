#!/usr/bin/env python

import os, sys

#print 'usage: %s [zp2mu_histos.root [histos.ps]' % sys.argv[0]    

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
    outputFile = 'histos.ps'

if not os.path.isfile(rootFile):
    sys.stderr.write('input file %s does not exist!\n' % rootFile)
    sys.exit(1)

print 'Creating %s from %s...' % (outputFile, rootFile)

f = TFile(rootFile)
histos = f.Zprime2muHistos
psd = PSDrawer(outputFile)

gStyle.SetHistMinimumZero(True)

def draw_page(page_type, histo_base_name, page_title, draw_log=False):
    if page_type == 'all':
        div = (2,5)
        levels = xrange(9)
    elif page_type == 'offline':
        div = (2,3)
        levels = xrange(4, 9)
    elif page_type == 'no_trig':
        div = (2,3)
        levels = [0] + range(4, 9)
    else:
        raise ValueError, 'page_type %s not recognized' % page_type
    
    pad = psd.new_page(page_title, div)
    for i, level in enumerate(levels):
        subpad = pad.cd(i+1)
        if draw_log:
            subpad.SetLogy(1)
        h = histos.Get('%s%i' % (histo_base_name, level))
        if h is not None:
            h.Draw()

def draw(name, opt=''):
    h = histos.Get(name)
    if h:
        h.Draw(opt)
    return h


pad = psd.new_page('Trigger information', (3,3))
for rec in xrange(1,4):
    offset = 3*(rec-1)
    pad.cd(offset+1).SetLogy(1)
    draw('TriggerBits%i' % rec)
    pad.cd(offset+2).SetLogy(1)
    draw('NLeptonsTriggered%i' % rec)
    pad.cd(offset+3).SetLogy(1)
    draw('NLeptonsFailed%i' % rec)

all_levels_pages = [
    ('NLeptons',      '# leptons/event', True),
    ('LeptonEta',     'Lepton #eta'),
    ('LeptonRap',     'Lepton rapidity'),
    ('LeptonPhi',     'Lepton #phi'),
    ('LeptonPt',      'Lepton pT'),
    ('LeptonPz',      'Lepton pz'),
    ('LeptonP',       'Lepton p'),
    ('LeptonPVsEta',  'Lepton p vs. #eta'),
    ('LeptonPtVsEta', 'Lepton pT vs. #eta'),
    ]

offline_pages = [
    ('IsoSumPt',      'Lepton isolation (#DeltaR < 0.3) #SigmapT', True),
    ('IsoEcal',       'Lepton isolation (#DeltaR < 0.3) ECAL', True),
    ('IsoHcal',       'Lepton isolation (#DeltaR < 0.3) HCAL', True),
    ('IsoNTracks',    'Lepton isolation (#DeltaR < 0.3) nTracks', True),
    ('IsoNJets',      'Lepton isolation (#DeltaR < 0.3) nJets', True),
    ('NPxHits',       'Lepton track # pixel hits', True),
    ('NStHits',       'Lepton track # strip hits', True),
    ('NTkHits',       'Lepton track # tracker hits', True),
    ('NMuHits',       'Lepton track # muon hits', True),
    ('NHits',         'Lepton track # total hits', True),
    ('Chi2dof',       'Lepton track #chi^{2}/dof'),
    ('TrackD0',       'Lepton track |d0|'),
    ('TrackDz',       'Lepton track |dz|')
    ]

no_trig_pages = [
    ('NDileptons',              '# dileptons/event', True),
    ('DileptonEta',             'Dilepton #eta'),
    ('DileptonRap',             'Dilepton rapidity'),
    ('DileptonPhi',             'Dilepton #phi'),
    ('DileptonPt',              'Dilepton pT'),
    ('DileptonPz',              'Dilepton pz'),
    ('DileptonP',               'Dilepton p'),
    ('DileptonPVsEta',          'Dilepton p vs. #eta'),
    ('DileptonPtVsEta',         'Dilepton pT vs. #eta'),
    ('DileptonMass',            'Dilepton mass'),
    ('DileptonWithPhotonsMass', 'Dilepton+photons mass'),
    ('DileptonSigns',           'Dilepton signs', True),
    ('DileptonDeltaPt',         'Dilepton |pT^{1}| - |pT^{2}|'),
    ('DileptonDeltaP',          'Dilepton |p^{1}| - |p^{2}|'),
    ('DileptonPtErrors',        'Dilepton #sigma_{pT}^{1} v. #sigma_{pT}^{2}')
    ]

for page_info in all_levels_pages:
    draw_page('all', *page_info)
for page_info in offline_pages:
    draw_page('offline', *page_info)
for page_info in no_trig_pages:
    draw_page('no_trig', *page_info)

psd.close()
f.Close()
