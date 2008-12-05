#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
    print 'usage: %s input.root [output.ps]' % sys.argv[0]    
    sys.exit(1)

# Run ROOT in batch mode.
sys.argv.append('-b')
from ROOT import *
sys.argv.remove('-b')

from PSDrawer import PSDrawer
    
rootFile = sys.argv[1]
if len(sys.argv) > 2:
    outputFile = sys.argv[2]
else:
    outputFile = 'output.ps'

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

all_levels_pages = [
    ('NLeptons',      '# leptons/event', True),
    ('LeptonEta',     'Lepton #eta'),
    ('LeptonRap',     'Lepton rapidity'),
    ('LeptonPhi',     'Lepton #phi'),
    ('LeptonPt',      'Lepton p_{T}'),
    ('LeptonPz',      'Lepton p_{z}'),
    ('LeptonP',       'Lepton p'),
    ('LeptonPVsEta',  'Lepton p vs. #eta'),
    ('LeptonPtVsEta', 'Lepton p_{T} vs. #eta'),
    ]

offline_pages = [
    ('IsoSumPt',      'Lepton isolation (#DeltaR < 0.3) #Sigmap_{T}', True),
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
    ('TrackD0',       'Lepton track |d_{0}|'),
    ('TrackDz',       'Lepton track |dz|')
    ]

no_trig_pages = [
    ('NDileptons',              '# dileptons/event', True),
    ('DileptonEta',             'Dilepton #eta'),
    ('DileptonRap',             'Dilepton rapidity'),
    ('DileptonPhi',             'Dilepton #phi'),
    ('DileptonPt',              'Dilepton p_{T}'),
    ('DileptonPz',              'Dilepton p_{z}'),
    ('DileptonP',               'Dilepton p'),
    ('DileptonPVsEta',          'Dilepton p vs. #eta'),
    ('DileptonPtVsEta',         'Dilepton p_{T} vs. #eta'),
    ('DileptonMass',            'Dilepton mass'),
    ('DileptonWithPhotonsMass', 'Dilepton+photons mass'),
    ('DileptonSigns',           'Dilepton signs', True),
    ('DileptonDeltaPt',         'Dilepton |p_{T}^{1}| - |p_{T}^{2}|'),
    ('DileptonDeltaP',          'Dilepton |p^{1}| - |p^{2}|'),
    ('DileptonPtErrors',        'Dilepton #sigma_{p_{T}}^{1} v. #sigma_{p_{T}}^{2}')
    ]

for page_info in all_levels_pages:
    draw_page('all', *page_info)
for page_info in offline_pages:
    draw_page('offline', *page_info)
for page_info in no_trig_pages:
    draw_page('no_trig', *page_info)

psd.close()
f.Close()
