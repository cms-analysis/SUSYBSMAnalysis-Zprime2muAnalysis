import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      #WJET
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0078400B-12FD-E411-8407-0025904C6374.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/00BB2319-15FD-E411-B35C-002590DBDFE0.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/00CE8D69-1FFD-E411-8737-0025905A6118.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/06CF92A0-F5FC-E411-9253-E0CB4E4408DD.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/06D587E7-03FD-E411-AA4B-001E67396905.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/082BF333-21FD-E411-95BB-002590DB9258.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0A7866E2-2DFD-E411-9E2E-002590DB9278.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0AD13C07-1BFD-E411-AB39-0025905A6118.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0E2733B5-E4FD-E411-BD3A-0CC47A0AD704.root',
#                                      '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/107B83F0-03FD-E411-BB01-000F53273500.root',
                                      #DY50
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0033A97B-8707-E511-9D3B-008CFA1980B8.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/003826EE-8807-E511-9628-008CFA06470C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/025EA7DE-9907-E511-B054-002590593876.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/06C07987-9A07-E511-B2C9-0025905A60D2.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/081564DD-9907-E511-B3EA-002590596498.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/08B5061B-9C07-E511-87F6-003048FFCB96.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0A691BA0-9207-E511-B3BA-008CFA064774.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0AFC520A-8607-E511-B099-001517FB20EC.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0C4979A5-9D07-E511-9551-0025905938D4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0C4DF708-7C07-E511-A290-00259073E3AE.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0E2E857A-8F07-E511-BE8B-001E675A6653.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0E6E5801-9B07-E511-8190-0025905A60EE.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0E849D70-8A07-E511-9B27-002590A887F2.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/1212C2E8-7E07-E511-A672-008CFA0A5844.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/121D820F-9C07-E511-854C-0025905A48E4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/125CA02D-9407-E511-8545-00259074AEDE.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/1465A073-8907-E511-B7D0-F45214C748D2.root',
                                      #DY100-200
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/18C4D462-7B2A-E511-A2D2-B499BAAC04E6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/22F6D18D-3A2A-E511-AF28-002354EF3BE0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/2AAB8D2B-792A-E511-AB9D-B499BAAC04F0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/300B0F2B-792A-E511-90ED-B499BAAC039C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/34E63DF8-0C2A-E511-84F4-A0040420FE80.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/5A624462-772A-E511-BF60-B499BAAC04E6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/92D13B26-792A-E511-A763-6C3BE5B5A038.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/C2213C5F-472A-E511-9A82-0002C92A1024.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/C2336FF2-0F2A-E511-9960-A0040420FE80.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/CC9C19E7-112A-E511-A7EF-B8CA3A70A5E8.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/E4A6C82B-742A-E511-99CB-B499BAAC04E6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/F60E86E9-792A-E511-A688-B499BAAC04E6.root'
                                      #DY200-400
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/101330DE-BF27-E511-BAEC-003048FFCB6A.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/1E686152-5728-E511-96F0-02163E011524.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/402930AB-D427-E511-AF3D-02163E01466D.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/54852789-C127-E511-8364-0025905A6122.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/7C837FA7-B627-E511-AA01-003048FFCC2C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/860A9276-5728-E511-8A5C-0025905A60EE.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/9EA68010-B727-E511-A5C5-0025905964B4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/A41DF1DD-BF27-E511-BBF9-0025905A60CE.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/B2C35A7D-5728-E511-A288-0025905A6094.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/C63F91EC-C227-E511-B334-0025905964C2.root'
                                      #DY400-500
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/0453FDA6-572A-E511-ACC8-002590A4FFE8.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/AC0341A6-572A-E511-AC29-002590200AF4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/EA7273A3-572A-E511-B128-001E67397756.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/EC3C2C04-FA29-E511-AC64-0002C90A3476.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/0C7F4A31-5C29-E511-9C52-0025B3E06378.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/1C9CE7C7-5F29-E511-BEDE-0002C94CD2C6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/247A4A2D-5C29-E511-B56C-002481E14F8C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/443B4EBA-902D-E511-A528-001E67396D10.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/4CD0D160-9728-E511-BE40-00259073E4CC.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/6094FC3C-9528-E511-B036-A0369F3102B6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/6E9B7B0A-D329-E511-93FC-20CF3027A633.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/9C2DD130-D329-E511-B1D2-0002C90F8030.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/B6F889E0-4E29-E511-893E-F45214CEF24A.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/CA545029-5C29-E511-BFA2-0025B3E05E1C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/D8270D61-4F29-E511-BC77-A0040420FE80.root'
                                      #DY500-700
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/1EA8037C-872D-E511-B22D-00266CFAE464.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/2AF3B8E3-792B-E511-8149-20CF305B0599.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/38F9B6F3-872D-E511-87E0-00266CF27170.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/4212435C-682B-E511-945E-0CC47A4DEE66.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/7AEAF41F-602B-E511-985E-0CC47A4D9A2E.root',
##                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/7C63DE4A-022C-E511-9992-002590A36084.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/80757942-822B-E511-8350-00259073E514.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/86CAE9E7-872D-E511-B7B2-3417EBE64C51.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/86D95D72-592B-E511-B22F-0CC47A4DEDB8.root',
##                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/9424E042-022C-E511-9D8E-002590A52B4A.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/B04413D8-622B-E511-9ED9-0CC47A4D9A2E.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/50000/B8A03121-602B-E511-8007-00259074AE98.root',
                                      #DY700-800
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/0E416950-FA29-E511-9A4B-549F358EB7CA.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/0EC928F0-D129-E511-B7B1-008CFA110AA8.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/3C3D3E3A-D529-E511-A73D-008CFA110AB4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/72BE72B3-D929-E511-907C-008CFA197A60.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/922C414F-CD29-E511-AA07-00259055220A.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/965800AB-DC29-E511-B229-008CFA197A60.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/96655751-EB29-E511-B3C9-008CFA197964.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/A82D35CD-A329-E511-A735-00266CFFA240.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/AAF0DB03-F729-E511-8A68-008CFA197410.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/BE3D3582-F029-E511-BD94-549F35AC7E8A.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/C6735DC9-E329-E511-9479-008CFA1979A0.root'
                                      #DY800-1000
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/0A0C6F31-F62A-E511-A9C3-0025B3E063AC.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/0EE82ED8-DA27-E511-836D-0025905A60E0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/129412B7-9F27-E511-BBA0-0025905B8576.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/1A880C81-A327-E511-BB49-0025905822B6.root',
##                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/902BDAE7-A027-E511-8E61-003048FFD7D4.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/A232159C-9827-E511-B8EE-0025905B85D6.root',
##                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/B61AA65A-A227-E511-81CA-0025905B85B2.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/C4DCE5EE-DA27-E511-80CB-0025905B855E.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/C82CC078-9627-E511-BDF5-0025905A60E0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/CA7D7549-9F27-E511-B0ED-003048FFD754.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/D4EF0DE0-DA27-E511-B74D-A0040420FE80.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/E0199D2E-F62A-E511-9597-002590200B00.root'
                                      ####DY1000-1500
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/08D109A3-602D-E511-BB4E-20CF3019DEFF.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/0AB70AFA-F129-E511-8E27-0CC47A4DE052.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/206791F6-EE29-E511-B206-0CC47A4DEEF0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/383614AE-602D-E511-A66B-00221982B698.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/8A23353C-D529-E511-89BA-002590775016.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/9636BE9A-412A-E511-8944-0CC47A4D99A0.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/96AA49A7-EA29-E511-B538-0CC47A4DEED6.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/983D9CF4-E729-E511-9788-20CF3027A5D5.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/A875CEA2-412A-E511-AEF8-00215AEDFD12.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/A8D15B32-612D-E511-8AB8-F04DA2770C8E.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/B05149E9-C429-E511-935C-00259073E44C.root',
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/30000/F0A8B49C-602D-E511-9133-20CF305B058D.root'
                                      #DY1500-2000
#                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/9A25A8BB-902D-E511-87FC-0025B3E06378.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/2CD1F837-7F2A-E511-BA03-000F532734B4.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/3CFE0B33-F72A-E511-A956-001E67396D51.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/3E521C7A-5E2A-E511-B9C9-000F532734A0.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/6E072944-612A-E511-963B-000F532734A0.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/70A89E35-F72A-E511-8FBE-002590A371AE.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/8E9007D5-F72A-E511-A2B2-002590A831CA.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/A22B6A53-632A-E511-B69A-B083FED76520.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/B0571BDF-642A-E511-87EF-AC853D9F5120.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/D68B814A-7F2A-E511-ABB7-B083FED76637.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/E6D524CB-792A-E511-B1FC-000F53273500.root',
                                      '/store/mc/RunIISpring15DR74/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/F416988B-832A-E511-8F6D-002590D9D9DA.root'
    )
)

process.maxEvents.input = -1


process.demo = cms.EDAnalyzer("WeightAnalyzer",

	#GenEventInfo = cms.InputTag("GenEventInfoProduct_generator","","SIM")
         genEventInfo = cms.InputTag('generator'),
        
	
)


process.p = cms.Path(process.demo)
