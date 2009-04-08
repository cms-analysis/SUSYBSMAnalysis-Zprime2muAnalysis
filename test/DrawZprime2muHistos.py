#!/usr/bin/env python

# usage: DrawZprime2muHistos.py [zp2mu_histos.root [histos.ps]

import os, sys
from SUSYBSMAnalysis.Zprime2muAnalysis.PSDrawer import PSDrawer
from ROOT import TFile
    
if len(sys.argv) > 1: rootFile = sys.argv[1]
else:                 rootFile = 'zp2mu_histos.root'

if len(sys.argv) > 2:
    outputFile = sys.argv[2]
else:
    outputFile = os.path.basename(rootFile).replace('.root', '.histos.ps')

if not os.path.isfile(rootFile):
    sys.stderr.write('input file %s does not exist!\n' % rootFile)
    sys.exit(1)

print 'Creating %s from %s...' % (outputFile, rootFile)

f = TFile(rootFile)
histos = f.Zprime2muHistos
psd = PSDrawer(outputFile)

################################################################################

pad = psd.new_page('Trigger information', (3,3))
for rec in xrange(psd.TRIG_START, psd.OFFLINE_START):
    offset = 3*(rec-1)
    pad.cd(offset+1).SetLogy(1)
    psd.draw_if(histos, 'TriggerBits%i' % rec)
    pad.cd(offset+2).SetLogy(1)
    psd.draw_if(histos, 'NLeptonsTriggered%i' % rec)
    pad.cd(offset+3).SetLogy(1)
    psd.draw_if(histos, 'NLeptonsFailed%i' % rec)

psd.rec_level_page(histos, 'all', 'NLeptons',      '# leptons/event', log_scale=True)
psd.rec_level_page(histos, 'all', 'LeptonEta',     'Lepton #eta')
#psd.rec_level_page(histos, 'all', 'LeptonRap',     'Lepton rapidity')
psd.rec_level_page(histos, 'all', 'LeptonPhi',     'Lepton #phi')
psd.rec_level_page(histos, 'all', 'LeptonPt',      'Lepton pT')
psd.rec_level_page(histos, 'all', 'LeptonPz',      'Lepton pz')
psd.rec_level_page(histos, 'all', 'LeptonP',       'Lepton p')
psd.rec_level_page(histos, 'all', 'LeptonPVsEta',  'Lepton p vs. #eta')
psd.rec_level_page(histos, 'all', 'LeptonPtVsEta', 'Lepton pT vs. #eta')

psd.rec_level_page(histos, 'offline', 'IsoSumPt',   'Lepton isolation (#DeltaR < 0.3) #SigmapT', log_scale=True)
psd.rec_level_page(histos, 'offline', 'IsoEcal',    'Lepton isolation (#DeltaR < 0.3) ECAL',     log_scale=True)
psd.rec_level_page(histos, 'offline', 'IsoHcal',    'Lepton isolation (#DeltaR < 0.3) HCAL',     log_scale=True)
psd.rec_level_page(histos, 'offline', 'IsoNTracks', 'Lepton isolation (#DeltaR < 0.3) nTracks',  log_scale=True)
psd.rec_level_page(histos, 'offline', 'IsoNJets',   'Lepton isolation (#DeltaR < 0.3) nJets',    log_scale=True)

psd.rec_level_page(histos, 'offline', 'NPxHits', 'Lepton track # pixel hits',   log_scale=True)
psd.rec_level_page(histos, 'offline', 'NStHits', 'Lepton track # strip hits',   log_scale=True)
psd.rec_level_page(histos, 'offline', 'NTkHits', 'Lepton track # tracker hits', log_scale=True)
psd.rec_level_page(histos, 'offline', 'NMuHits', 'Lepton track # muon hits',    log_scale=True)
psd.rec_level_page(histos, 'offline', 'NHits',   'Lepton track # total hits',   log_scale=True)

psd.rec_level_page(histos, 'offline', 'Chi2dof', 'Lepton track #chi^{2}/dof')
psd.rec_level_page(histos, 'offline', 'TrackD0', 'Lepton track |d0|')
psd.rec_level_page(histos, 'offline', 'TrackDz', 'Lepton track |dz|')

psd.rec_level_page(histos, 'no_trig', 'NDileptons',      '# dileptons/event', log_scale=True)
psd.rec_level_page(histos, 'no_trig', 'DileptonEta',     'Dilepton #eta')
psd.rec_level_page(histos, 'no_trig', 'DileptonRap',     'Dilepton rapidity')
psd.rec_level_page(histos, 'no_trig', 'DileptonPhi',     'Dilepton #phi')
psd.rec_level_page(histos, 'no_trig', 'DileptonPt',      'Dilepton pT')
psd.rec_level_page(histos, 'no_trig', 'DileptonPz',      'Dilepton pz')
psd.rec_level_page(histos, 'no_trig', 'DileptonP',       'Dilepton p')
psd.rec_level_page(histos, 'no_trig', 'DileptonPVsEta',  'Dilepton p vs. #eta')
psd.rec_level_page(histos, 'no_trig', 'DileptonPtVsEta', 'Dilepton pT vs. #eta')

psd.rec_level_page(histos, 'no_trig', 'DileptonMass',            'Dilepton mass')
psd.rec_level_page(histos, 'no_trig', 'DileptonWithPhotonsMass', 'Dilepton+photons mass')

psd.rec_level_page(histos, 'no_trig', 'DileptonSigns',    'Dilepton signs', log_scale=True)
psd.rec_level_page(histos, 'no_trig', 'DileptonDeltaPt',  'Dilepton |pT^{1}| - |pT^{2}|')
psd.rec_level_page(histos, 'no_trig', 'DileptonDeltaP',   'Dilepton |p^{1}| - |p^{2}|')
psd.rec_level_page(histos, 'no_trig', 'DileptonPtErrors', 'Dilepton #sigma_{pT}^{1} v. #sigma_{pT}^{2}')

psd.rec_level_page(histos, 'cocktail', 'TeVMuonCocktailSource', 'Original rec level of cocktail muons', draw_opt='HIST TEXT0', hist_cmds=[('SetMarkerSize', 1), ('SetStats', 0)])

psd.close()
f.Close()
