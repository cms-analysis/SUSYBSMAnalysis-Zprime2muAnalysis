import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysisCommon_cff import *
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muMassReach_cfi import *

files = [
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_1.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_2.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_3.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_4.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_5.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_6.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_7.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_8.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_9.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_10.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_11.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_12.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_13.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_14.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_15.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_16.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_17.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_18.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_19.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_20.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_21.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_22.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_23.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_24.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/ZPSSMmumu_M1200_Mcut600-MC_31X_V3/PYTHIA6_ZPSSMmumu_M1200_Mcut600_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_25.root',
#    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/Zssm1000_fm_RECO.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_1.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_2.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_3.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_4.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_5.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_6.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_7.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_8.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_9.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_10.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_11.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_12.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_13.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_14.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_15.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_16.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_17.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_18.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_19.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_20.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_21.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_22.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_23.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_24.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut500-MC_31X_V3/PYTHIA6_DYmumu_Mcut500_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_25.root',
#    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_1_6_7/dy_above200_RECO.root'
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_1.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_2.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_3.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_4.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_5.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_6.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_7.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_8.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_9.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_10.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_11.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_12.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_13.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_14.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_15.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_16.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_17.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_18.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_19.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_20.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_21.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_22.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_23.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_24.root',
    'rfio:/castor/cern.ch/user/t/tucker/CMSSW_3_1_2/DYmumu_Mcut200-MC_31X_V3/PYTHIA6_DYmumu_Mcut200_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_HLT8E29_RECO_MC_31X_V3_25.root'
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root',
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root',
#   'file:/data0/slava/data/ZPSSMmumu_M1000_Mcut400_1.root'
    ]

# process = makeZprime2muAnalysisProcess(files)
process = makeZprime2muAnalysisProcess(files, maxEvents=-1, skipPAT=True, disableElectrons=True, conditionsGlobalTag='STARTUP31X_V2::All') #, useTrigger=False)

#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/dy_for1TeV_official.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm1000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm1500_private_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1500_private_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi1500_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zssm2000_massreach.cfi")
#process.include("SUSYBSMAnalysis/Zprime2muAnalysis/test/getSamples/Zpsi2000_massreach.cfi")
    
attachMassReach(process)
process.Zprime2muMassReach.verbosity = 0

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


process.Zprime2muMassReach.dataSet       = "Zssm1200_EPE"

#process.Zprime2muMassReach.dataSet       = "Zpsi1000"
#process.Zprime2muMassReach.dataSet       = "Zssm1500"
#process.Zprime2muMassReach.dataSet       = "Zpsi1500"
#process.Zprime2muMassReach.dataSet       = "Zssm2000"
#process.Zprime2muMassReach.dataSet       = "Zpsi2000"
#process.Zprime2muMassReach.DYEvents      = True
#process.Zprime2muMassReach.BackgroundFit = True
#process.Zprime2muMassReach.FixedMass     = False
#process.Zprime2muMassReach.FixedFWHM     = False
#process.Zprime2muMassReach.ExpPlots      = True
#process.Zprime2muMassReach.BinnedFit     = True
process.Zprime2muMassReach.intLumi       = 0.05

