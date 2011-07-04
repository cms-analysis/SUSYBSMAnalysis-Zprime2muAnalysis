from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/PromptReco-v4/000/167/913/F2B6FB32-7AA3-E011-BCAF-BCAEC5329710.root']
process.GlobalTag.globaltag = 'GR_R_42_V13::All'
process.maxEvents.input = -1
process.options.wantSummary = True

process.CheckPrescale = cms.EDAnalyzer('CheckPrescale',
                                       hlt_process_name = cms.string('HLT'),
                                       trigger_paths = cms.vstring('HLT_Mu30_v1', 'HLT_Mu30_v2', 'HLT_Mu30_v3', 'HLT_Mu30_v4', 'HLT_Mu30_v5')
                                       )
process.p = cms.Path(process.CheckPrescale)
