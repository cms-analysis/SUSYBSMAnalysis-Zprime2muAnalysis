import os, FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.p = cms.Path(process.countPatLeptons)

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
#scheduler = %(scheduler)s
scheduler = remoteGlidein

use_server = 0

[CMSSW]
datasetpath = %(dataset)s
pset = %(pset)s
get_edm_output = 1
%(job_control)s

use_dbs3=1

[USER]
eMail = raffaella.radogna@cern.ch

ui_working_dir = crab/crab_datamc_%(name)s
copy_data = 1
#storage_element = T2_IT_Bari
storage_element = T2_IT_Legnaro
#storage_path = /lustre/cms/store/user/rradogna/
check_user_remote_dir = 0
publish_data = 1
publish_data_name = datamc_%(name)s
#dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet

return_data = 0
user_remote_dir = Zprime/crab_datamc_%(name)s

[GRID]
#ce_white_list = T2_EE_Estonia
'''

os.system('mkdir -p crab/psets')
