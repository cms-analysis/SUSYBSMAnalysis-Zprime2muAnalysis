#!/usr/bin/env python

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection_cff import loose_cut, trigger_match, tight_cut, allDimuons

# Define the numerators and denominators, removing cuts from the
# allDimuons maker. "NoX" means remove cut X (i.e. the denominators),
# "NoNo" means remove nothing (i.e. the numerators). This will break
# if loose_, tight_ cut strings are changed...
process.allDimuonsNoIso     = allDimuons.clone(loose_cut = loose_cut.replace(' && isolationR03.sumPt < 10', ''))
process.allDimuonsNoTkHits  = allDimuons.clone(loose_cut = loose_cut.replace(' && innerTrack.hitPattern.numberOfValidTrackerHits >= 10', ''))
process.allDimuonsNoDB      = allDimuons.clone(tight_cut = tight_cut.replace('dB < 0.2 && ', ''))
process.allDimuonsNoGlbChi2 = allDimuons.clone(tight_cut = tight_cut.replace(' && globalTrack.normalizedChi2 < 10', ''))
process.allDimuonsNoPxHits  = allDimuons.clone(tight_cut = tight_cut.replace(' && innerTrack.hitPattern.numberOfValidPixelHits >= 1', ''))
process.allDimuonsNoMuStns  = allDimuons.clone(tight_cut = tight_cut.replace(' && globalTrack.hitPattern.muonStationsWithValidHits >= 2', ''))
process.allDimuonsNoTkMuon  = allDimuons.clone(tight_cut = tight_cut.replace(' && isTrackerMuon', ''))
process.allDimuonsNoTrgMtch = allDimuons.clone(tight_cut = tight_cut.replace(' && ' + vbtf_trigger_match, ''))
process.allDimuonsNoNo      = allDimuons.clone()

alldimus = [x for x in dir(process) if 'allDimuonsNo' in x]
process.p = cms.Path(process.goodDataFilter * process.muonPhotonMatch * process.leptons * reduce(lambda x,y: x*y, [getattr(process, x) for x in alldimus]))

# For all the allDimuons producers, make dimuons producers, and
# analyzers to make the histograms.
for alld in alldimus:
    dimu = process.dimuons.clone(src = alld)
    name = alld.replace('allD', 'd')
    setattr(process, name, dimu)
    hists = HistosFromPAT.clone(dilepton_src = name, leptonsFromDileptons = True)
    setattr(process, name.replace('dimuons', ''), hists)
    process.p *= dimu * hists

# Handle the cuts that have to be applied at the
# Zprime2muCompositeCandidatePicker level.
process.p *= process.allDimuons
process.dimuonsNoB2B     = process.dimuons.clone(back_to_back_cos_angle_min = -1)
process.dimuonsNoVtxProb = process.dimuons.clone(vertex_chi2_max = -1)

for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb']:
    hists = HistosFromPAT.clone(dilepton_src = dimu, leptonsFromDileptons = True)
    setattr(process, dimu.replace('dimuons', ''), hists)
    process.p *= getattr(process, dimu) * hists


if 'olddata' in sys.argv:
    process.source.fileNames = ['file:work/daata/jul15.root', 'file:work/daata/prompt.root']
    process.TFileService.fileName = 'ana_nminus1_olddata.root'
elif 'data' in sys.argv:
    process.source.fileNames = ['file:../DataMCSpectraComparison/crab/crab_datamc_promptB_all/res/merged.root']
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*open('../DataMCSpectraComparison/ana_datamc/ana_datamc_data_promptB_allgood.cmssw').read().split(','))
    process.TFileService.fileName = 'ana_nminus1_data.root'
