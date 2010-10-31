# The below is copy-pasted from the DBS config for
# /DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO 
# and modified to do what we want.

import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('/store/user/tucker/dy120-HLT-384p3-START38_V12/dy120-HLT-384p3-START38_V12/9084f19c4f7cc164e43ee71169c6036c/hlt_1_1_7V0.root'))

process.RECOSIMoutput = cms.OutputModule('PoolOutputModule',
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('reco.root'),
    dataset = cms.untracked.PSet(filterName = cms.untracked.string(''), dataTier = cms.untracked.string('GEN-SIM-RECO'))
)

process.mix.playback = True
process.GlobalTag.globaltag = 'START38_V12::All'

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.RECOSIMoutput_step)

def customiseCommon(process):
    
    #####################################################################################################
    ####
    ####  Top level replaces for handling strange scenarios of early collisions
    ####

    ## TRACKING:
    ## Skip events with HV off
    process.newSeedFromTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
    process.newSeedFromPairs.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
    process.secTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
    process.fifthSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 20000
    process.fourthPLSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
    process.thTripletsA.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000
    process.thTripletsB.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000

    ###### FIXES TRIPLETS FOR LARGE BS DISPLACEMENT ######

    ### prevent bias in pixel vertex
    process.pixelVertices.useBeamConstraint = False
    
    ### pixelTracks
    #---- new parameters ----
    process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = 4.06
    process.pixelTracks.RegionFactoryPSet.RegionPSet.originHalfLength = cms.double(40.6)

    ### 0th step of iterative tracking
    #---- new parameters ----
    process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(4.06)  
    process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 40.6

    ### 2nd step of iterative tracking
    #---- new parameters ----
    process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)  
    process.secTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 44.7

    ## ECAL 
    process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]

    ###
    ###  end of top level replacements
    ###
    ###############################################################################################

    return (process)

process = customiseCommon(process)


import sys, os
if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = condor

[CMSSW]
datasetpath = %(dataset)s
dbs_url = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
pset = psets/reco_crab.py
get_edm_output = 1
total_number_of_events = -1
events_per_job = 5000

[USER]
ui_working_dir = crab/crab_reco_%(name)s
copy_data = 1
storage_element = T3_US_FNALLPC
check_user_remote_dir = 0
publish_data = 1
publish_data_name = %(name)s-RECO-384p3-START38_V12
dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
'''

    x = [
        '/dy120-HLT-384p3-START38_V12/tucker-dy120-HLT-384p3-START38_V12-9084f19c4f7cc164e43ee71169c6036c/USER',
        '/dy200-HLT-384p3-START38_V12/tucker-dy200-HLT-384p3-START38_V12-5f9cf42c6642247cfe6d8a7f720fef92/USER',
        '/dy500-HLT-384p3-START38_V12/tucker-dy500-HLT-384p3-START38_V12-66bfe43247a0156f0dd9742592a27e86/USER',
        '/dy800-HLT-384p3-START38_V12/tucker-dy800-HLT-384p3-START38_V12-39f7dba0ea9e98694efe81e32cc2b888/USER',
        '/zssm1000-HLT-384p3-START38_V12/tucker-zssm1000-HLT-384p3-START38_V12-358cd7657aecca27bbbd98c24f87a692/USER',
        '/zssm1250-HLT-384p3-START38_V12/tucker-zssm1250-HLT-384p3-START38_V12-3ddc24773901216a765835604ed443a9/USER',
        '/zssm1500-HLT-384p3-START38_V12/tucker-zssm1500-HLT-384p3-START38_V12-878ab2530af7e065fdb4b3bd86bfd641/USER',
        '/zssm1750-HLT-384p3-START38_V12/tucker-zssm1750-HLT-384p3-START38_V12-aabf1a021743d24ee2be2b03e685b061/USER',
        ]
    [
        '/zssm2000-HLT-384p3-START38_V12/tucker-zssm2000-HLT-384p3-START38_V12-a142f9a0730ff4a0ece0671b086562aa/USER',
        '/dy1200-HLT-384p3-START38_V12/tucker-dy1200-HLT-384p3-START38_V12-04d5b4bf88ca36f673706f20225f27e1/USER',
        '/dy1500-HLT-384p3-START38_V12/tucker-dy1500-HLT-384p3-START38_V12-8bc32f8a2989ed61ce7671e38b2d3fad/USER',
        '/dy1800-HLT-384p3-START38_V12/tucker-dy1800-HLT-384p3-START38_V12-cb77c0b470476268f62c726ea8aa1569/USER',
        '/zssm250-HLT-384p3-START38_V12/tucker-zssm250-HLT-384p3-START38_V12-75eefaca77eb8f49415074195ddb8263/USER',
        '/zssm500-HLT-384p3-START38_V12/tucker-zssm500-HLT-384p3-START38_V12-cbe4fb75afa4c5e1ffd3731ccf796faf/USER',
        '/zssm750-HLT-384p3-START38_V12/tucker-zssm750-HLT-384p3-START38_V12-bcc6ebe5127572de453dfa0a7c47305e/USER',
        '/zpsi250-HLT-384p3-START38_V12/tucker-zpsi250-HLT-384p3-START38_V12-1f7ca1ac17da2c994083ecf08488b44e/USER',
        '/zpsi500-HLT-384p3-START38_V12/tucker-zpsi500-HLT-384p3-START38_V12-1edd0bbd09aac0c894f11084b0ec51ba/USER',
        '/zpsi750-HLT-384p3-START38_V12/tucker-zpsi750-HLT-384p3-START38_V12-90bcb360b0ce3002494c18aceeefe6fa/USER',
        '/zpsi1000-HLT-384p3-START38_V12/tucker-zpsi1000-HLT-384p3-START38_V12-7e8542ec21a18ac30f61577cbb0b375e/USER',
        '/zpsi1250-HLT-384p3-START38_V12/tucker-zpsi1250-HLT-384p3-START38_V12-c87d32e77274f8c2d711d1b2f28852da/USER',
        '/zpsi1500-HLT-384p3-START38_V12/tucker-zpsi1500-HLT-384p3-START38_V12-f08d511d159435e69a92913e36a7a778/USER',
        '/zpsi1750-HLT-384p3-START38_V12/tucker-zpsi1750-HLT-384p3-START38_V12-684319d877b85045c3e5b994dd784674/USER',
        '/zpsi2000-HLT-384p3-START38_V12/tucker-zpsi2000-HLT-384p3-START38_V12-115dfa9e5b2cce4de155acb865eeabf2/USER',
        ]
        
    os.system('mkdir -p psets')
    os.system('mkdir -p crab')
    os.system('cp reco.py psets/reco_crab.py')

    for dataset in x:
        name = dataset.split('-')[0][1:]
        print name
        open('crab.cfg', 'wt').write(crab_cfg % locals())
        if not 'testing' in sys.argv:
            os.system('crab -create -submit')
        
        
