#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.ResolutionUsingMC_cfi import ResolutionUsingMC

levels = ['Global', 'TkOnly', 'TPFMS', 'Picky', 'TMR', 'TuneP', 'SigmaSwitch']

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import rec_levels
rec_levels(process, levels)

useMC = True
process.p = cms.Path(process.Zprime2muAnalysisSequence)
for level in levels:
    lepton_src = cms.InputTag('leptons', level)
    dilepton_src = cms.InputTag('dimuons' + level)
    
    h = HistosFromPAT.clone(lepton_src=lepton_src, dilepton_src=dilepton_src)
    setattr(process, 'HistosFromPAT' + level, h)

    if useMC:
        r = ResolutionUsingMC.clone(lepton_src=lepton_src, dilepton_src=dilepton_src)
        setattr(process, 'ResolutionUsingMC' + level, r)

    process.p *= h*r
