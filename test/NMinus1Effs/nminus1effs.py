#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff import vbtf_loose, vbtf_tight, allDimuons

# Define the numerators and denominators, removing cuts from the VBTF
# allDimuons maker. "NoX" means remove cut X (i.e. the denominators),
# "NoNo" means remove nothing (i.e. the numerators).
#
# This will break if vbtf_loose, vbtf_tight cut strings are changed...
#
process.allDimuonsVBTFNoEta24   = allDimuons.clone(loose_cut = vbtf_loose.replace(' && abs(eta) < 2.4', ''))
process.allDimuonsVBTFNoIso3    = allDimuons.clone(loose_cut = vbtf_loose.replace(' && isolationR03.sumPt < 3', ''))
process.allDimuonsVBTFNoTkHits  = allDimuons.clone(loose_cut = vbtf_loose.replace(' && innerTrack.hitPattern.numberOfValidTrackerHits >= 10', ''))
process.allDimuonsVBTFNoDB      = allDimuons.clone(tight_cut = vbtf_tight.replace('dB < 0.2 && ', ''))
process.allDimuonsVBTFNoGlbChi2 = allDimuons.clone(tight_cut = vbtf_tight.replace(' && globalTrack.normalizedChi2 < 10', ''))
process.allDimuonsVBTFNoPxHits  = allDimuons.clone(tight_cut = vbtf_tight.replace(' && innerTrack.hitPattern.numberOfValidPixelHits >= 1', ''))
process.allDimuonsVBTFNoMuStns  = allDimuons.clone(tight_cut = vbtf_tight.replace(' && globalTrack.hitPattern.muonStationsWithValidHits >= 2', ''))
process.allDimuonsVBTFNoTkMuon  = allDimuons.clone(tight_cut = vbtf_tight.replace(' && isTrackerMuon', ''))
process.allDimuonsVBTFNoEta21   = allDimuons.clone(tight_cut = vbtf_tight.replace(' && abs(eta) < 2.1', ''))
process.allDimuonsVBTFNoTrgMtch = allDimuons.clone(tight_cut = vbtf_tight.replace(' && !triggerObjectMatchesByPath("HLT_Mu9").empty() && abs(triggerObjectMatchesByPath("HLT_Mu9").at(0).eta()) < 2.1', ''))
process.allDimuonsVBTFNoIso10   = process.allDimuonsVBTFNoIso3.clone()
process.allDimuonsVBTFNoNo      = allDimuons.clone()
process.allDimuonsVBTFNoNoIso10 = allDimuons.clone(loose_cut = vbtf_loose.replace(' && isolationR03.sumPt < 3', ' && isolationR03.sumPt < 10'))

alldimus = [x for x in dir(process) if 'allDimuonsVBTFNo' in x]
process.p = cms.Path(process.goodDataFilter * process.hltFilter * process.muonPhotonMatch * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))

# For all the allDimuons producers, make dimuons producers, and
# analyzers to make the histograms.
for alld in alldimus:
    dimu = process.dimuons.clone(src = alld)
    name = alld.replace('allD', 'd')
    setattr(process, name, dimu)
    hists = HistosFromPAT.clone(dilepton_src = name, leptonsFromDileptons = True)
    setattr(process, name.replace('dimuons', ''), hists)
    process.p *= dimu * hists

################################################################################

from SUSYBSMAnalysis.Zprime2muAnalysis.tools import files_from_dbs
if 'data' in sys.argv:
    process.source.fileNames = ['file:work/daata/jul15.root', 'file:work/daata/prompt.root']
    process.TFileService.fileName = 'ana_nminus1_data.root'
elif 'zp1000' in sys.argv:
    process.source.fileNames = files_from_dbs('/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/tucker-dyzpforeff_zp1000-9caa3d7638ff33984d7b458a78e3e8dd/USER')
    process.p.remove(process.goodDataFilter)
    process.TFileService.fileName = 'ana_nminus1_zp1000.root'
elif 'dy120' in sys.argv:
    process.source.fileNames = files_from_dbs('/DYToMuMu_M-120_7TeV-pythia6/tucker-dyzpforeff_dy120-9caa3d7638ff33984d7b458a78e3e8dd/USER')
    process.p.remove(process.goodDataFilter)
    process.TFileService.fileName = 'ana_nminus1_dy120.root'
elif 'zmumu' in sys.argv:
    process.source.fileNames = files_from_dbs('/Zmumu_M20_CTEQ66-powheg/tucker-datamc_zmumu-ea1b4401edd0c9e8af9e80917519ee4e/USER')
    process.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'REDIGI36X')
    process.TFileService.fileName = 'ana_nminus1_zmumu.root'
elif 'ttbar' in sys.argv:
    process.source.fileNames = files_from_dbs('/TTbarJets-madgraph/tucker-datamc_ttbar-f1606b2f2aa07e70082d78e786896133/USER')
    process.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'REDIGI')
    process.TFileService.fileName = 'ana_nminus1_ttbar.root'
