import os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.p = cms.Path(process.type0PFMEtCorrection * process.patDefaultSequence)

# Loose cut on muons; stronger cuts to be applied for different
# sets of plots (e.g. add our isolation cut, or apply VBTF).
process.selectedPatMuons.cut = 'isGlobalMuon && pt > 20'

# Want to select only events that have at least two leptons (=
# muons+electrons), where the electrons must pass HEEP id, but don't
# want to force HEEP id on selectedPatElectrons so as not to screw up
# the jet cleaning until we study this.
process.heepPatElectrons = cms.EDFilter('PATElectronSelector',
                                        src = cms.InputTag('patElectrons'),
                                        cut = cms.string('userInt("HEEPId") == 0')
                                        )
process.patDefaultSequence.replace(process.selectedPatElectrons, process.selectedPatElectrons * process.heepPatElectrons)
process.countPatMuons.minNumber = 0
process.countPatLeptons.electronSource = cms.InputTag('heepPatElectrons')
process.countPatLeptons.minNumber = 2

crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = %(scheduler)s

[CMSSW]
datasetpath = %(dataset)s
pset = %(pset)s
get_edm_output = 1
%(job_control)s

[USER]
ui_working_dir = crab/crab_datamc_%(name)s
copy_data = 1
storage_element = T3_US_FNALLPC
check_user_remote_dir = 0
publish_data = 1
publish_data_name = datamc_%(name)s
dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet

[GRID]
#ce_white_list = T2_EE_Estonia
'''

os.system('mkdir -p crab/psets')
