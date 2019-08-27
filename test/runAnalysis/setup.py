#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =['dummyFile']

process.maxEvents.input = -1
isMC = %(isMC)s
addNTuples = %(addNTuples)s
year = %(year)d
sampleName = '%(name)s'
process.GlobalTag.globaltag = '%(GT)s'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000
#import for high pT muon triggers
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, trigger_paths, overall_prescale, offline_pt_threshold, trigger_filters, trigger_filters2016, trigger_path_names, trigger_path_names2016, trigger_match_2018, trigger_match_2016,trigger_path_full_names, trigger_path_full_names2016
#import for prescaled low pT muon triggers
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import prescaled_trigger_pt_threshold, overall_prescale_2016, overall_prescale_2017, overall_prescale_2018, prescaled_trigger_filters_16, prescaled_trigger_path_names_16, prescaled_trigger_path_full_names_16, prescaled_trigger_match_2016, prescaled_trigger_filters_18, prescaled_trigger_path_names_18, prescaled_trigger_path_full_names_18, prescaled_trigger_match_2018, prescaled_trigger_path_name_list_16, prescaled_trigger_path_name_list_17, prescaled_trigger_path_name_list_18, prescaled_offline_pt_threshold

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)
