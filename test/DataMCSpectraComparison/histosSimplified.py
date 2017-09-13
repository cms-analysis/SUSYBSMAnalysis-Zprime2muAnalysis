#!/usr/bin/env python


Electrons = False

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =[#'file:./pat.root'
'/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/2810A5DC-03C8-E611-B20C-001E67504B25.root',
# '/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/824C363B-0AC8-E611-B4A5-20CF3027A580.root',
# '/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/243D09B4-90D1-E611-B0FA-001E674DA347.root',
# '/store/mc/RunIISummer16MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06BD3929-EDC7-E611-8EC3-02163E019C96.root',

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/0C7F3206-4941-E711-B557-0242AC130003.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/46D57269-1442-E711-B0AB-F4E9D497BBE0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/5C0BBCA6-1442-E711-BC65-0025907DE22C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/8084889C-1442-E711-9BA9-0242AC130002.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/84525990-1442-E711-A6F6-0025901C1876.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/8E688071-1442-E711-AF3A-FA163E3BDB32.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/988C4BF1-EB41-E711-BF9A-001E67DFF609.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/98A57765-1442-E711-AA97-FA163E49A90A.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/A677C09B-1442-E711-9BB4-B083FED00118.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/CC28BA7E-E341-E711-BCDC-B8CA3A709028.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/F867027C-1442-E711-8EF3-0CC47A78A3B4.root",

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/62DF9A66-1B3B-E711-AFFB-3417EBE5291B.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/82607D92-1B3B-E711-92E4-00266CFFCB28.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/9445D9C2-9E39-E711-AAB6-003048CB8652.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/DC1E16AA-1B3B-E711-BD61-68B59972BFD8.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/164D5782-CF37-E711-86A8-44A842CF058B.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/2EB31D0C-CF37-E711-8F04-FA163E148722.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/4C2C7C07-CF37-E711-9A6F-34E6D7E387A8.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/4EF1D9B4-C537-E711-90F2-008CFA197A60.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/604AED09-CF37-E711-990A-0CC47A745250.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/7A70B620-AA37-E711-A64C-3417EBE2F1CF.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/82AF03A0-9637-E711-BE4B-A0369F7FE960.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/82D75209-CF37-E711-9766-00259021A262.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/8A970B0B-CF37-E711-BF24-D4AE526A0419.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/B4719905-CF37-E711-B007-848F69FF914E.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/B825F4FC-CE37-E711-BC83-0025907277CE.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/C051E405-CF37-E711-830B-0242AC130002.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/E2A1F1D7-E736-E711-9E33-14187741120B.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/EA7BF004-CF37-E711-8126-FA163E57BB39.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/F07C4503-CF37-E711-9BBC-A0369F83627E.root",

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/2C3F9634-D736-E711-B084-002590D9D8BA.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/40639EE5-D636-E711-A893-24BE05CEEB81.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/561676E5-D636-E711-8F9D-549F3525DD6C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/72D5CFD3-D636-E711-8D97-0025905B8574.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/86F0AACE-D636-E711-9F19-3417EBE64BE5.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/8E94E9DE-D636-E711-8895-0242AC110004.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/9E3866CC-D636-E711-BA5B-0025905C43EA.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/CC6788CF-D636-E711-9865-FA163E726FD0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/D60EC4D3-D636-E711-85E6-00266CFFBDB4.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/E81A7A2C-D736-E711-86A1-0242AC130002.root",

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1250FF54-CF36-E711-87AD-0025905C3E36.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/14012D6A-CF36-E711-A4C4-0242AC110004.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A882B57-CF36-E711-AC33-0CC47A7C349C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1E235B5B-CF36-E711-A8E1-C4346BC7EDD8.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/32A5F25A-CF36-E711-8C4D-A0369F7FC934.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/562FDC68-CF36-E711-B710-7CD30ACE176E.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/602EA087-CF36-E711-A568-549F3525A64C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/9A092B6A-CF36-E711-A281-001E6739689C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/9C188F54-CF36-E711-950B-24BE05C44BB1.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/A6017B93-CF36-E711-9BD0-00259075D702.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/B4A8D85D-CF36-E711-B479-0CC47A4D7654.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/BAE42461-CF36-E711-A592-002590FD5E80.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/EC55B974-CF36-E711-8CFA-002590D60062.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/EE98F958-CF36-E711-80A6-A4BF0102572F.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/F0BCCF48-CF36-E711-8E40-FA163E5A1C0E.root",

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/007EBC7D-8C36-E711-B6AC-44A842CFD65A.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/16E8177F-8C36-E711-9879-0CC47A0107D0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/3201DD2F-8C36-E711-A35E-C4346BC76CD8.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/4C950A2C-8C36-E711-B07F-0CC47A4D76D2.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/5CF0C97E-8C36-E711-A931-C0BFC0E5685E.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/5E38007C-8C36-E711-9FED-008CFAF72A64.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/9CC84185-8C36-E711-ADAB-0026B92779FE.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/B820D181-8C36-E711-B603-FA163ED8B1C5.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/C2607881-8C36-E711-BB8B-901B0E6459E0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/CE51747E-8C36-E711-8A32-A0369FC5EB28.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/DE221930-8C36-E711-B115-48D539F3867C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/DE98887C-8C36-E711-8DCA-B083FED13803.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/ECC19231-8C36-E711-9739-3417EBE722FA.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/FC247932-8C36-E711-BF1C-0025905B85AE.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/FC3CBC3B-8C36-E711-9529-0242AC11000B.root",

#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/062F8710-D236-E711-A706-14DDA90900BB.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/069864E8-D136-E711-ADFC-008CFA197AA0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/148852D9-D136-E711-A77E-B083FED177B1.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/1A115C09-D236-E711-99E6-008CFAFBF2CA.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/30C050E6-D136-E711-92CE-485B3919F111.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/389D33CE-D136-E711-9095-008CFA1C907C.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/489242BF-D136-E711-93AD-001E677926C0.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/500591CD-D136-E711-8106-24BE05CE2EC1.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/96512B0F-D236-E711-B999-20CF3027A57D.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/B69A79C6-D136-E711-B1A5-0025905A60B2.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/BCFB8110-D236-E711-9051-B083FED13803.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/C4C8F5BD-D136-E711-BA56-002590D0AF4A.root",
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/EE52C4CC-D136-E711-9695-0242AC110002.root,"
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/F6151DE0-D136-E711-B940-001E67457DFA.root,"
#"/store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVDesRR_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/F66F7DC1-D136-E711-A31F-7CD30ACE17F2.root,"


#            '/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/162AD1DB-1E98-E611-9893-008CFA56D58C.root',
#            '/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/1A1F07FF-2698-E611-915C-0242AC130004.root',
                           

			   
			   ]
process.maxEvents.input = -1
isMC = False
for fileName in process.source.fileNames:
	if "Run2016H" in fileName:
		process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14'
	elif "Run2016" in fileName:
		process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'
	else:
		process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
		isMC = True
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)

from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor on DY. If used, the k factor will be applied to ALL samples submitted. #####

	
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
#import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2011EPS_cff as OurSelection2011EPS
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016






# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
	'Our2016'  : OurSelection2016,
	'Simple'   : OurSelection2016, # The selection cuts in the module listed here are ignored below.
	}

if Electrons:
	dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
		('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
		('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
		('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
		('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
		('MuonsPlusElectronsMinus',      '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == -2'),
		('MuonsMinusElectronsPlus',      '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == 2'),
		('MuonsPlusElectronsPlus',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == -24'),
		('MuonsMinusElectronsMinus',     '%(leptons_name)s:muons@- %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == 24'),
		('MuonsElectronsOppSign',        '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     ''),
		('MuonsElectronsSameSign',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
		('MuonsElectronsAllSigns',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
		]
	
	cuts = {
		'Our2016'  : OurSelection2016,
		'EmuVeto'  : OurSelection2016, # this switches on the dRMuEl veto
		'Simple'   : OurSelection2016, # The selection cuts in the module listed here are ignored below.
		}
	

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
	# Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present
    path_list.append(process.egmGsfElectronIDSequence)
	    
   
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut
    leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)

    if Electrons:
	    if cut_name == 'EmuVeto':
		    leptons.electron_muon_veto_dR = 0.1

    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'MuPrescaled' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
            assert alldil.tight_cut == trigger_match
            alldil.tight_cut = prescaled_trigger_match

    # Histos now just needs to know which leptons and dileptons to use.
      
	histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)


       #define the list of MC samples to be read here. be careful that if WWinclusive or tautau sample are not commented it will apply the filters when running locally.

    samples = [
	('dy50to120', '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
    	('dy120to200', '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
	('dy200to400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
    	('dy400to800', '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
	('dy800to1400', '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
    	('dy1400to2300', '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
	('dy2300to3500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
    	('dy3500to4500', '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
	('dy4500to6000', '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#
     	('WZ', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('WZ_ext', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
#     	
 	('ZZ', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('ZZ_ext', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
# 
     	('WWinclusive', '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('WW200to600', '/WWTo2L2Nu_Mll_200To600_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('WW600to1200', '/WWTo2L2Nu_Mll_600To1200_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('WW1200to2500', '/WWTo2L2Nu_Mll_1200To2500_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('WW2500', '/WWTo2L2Nu_Mll_2500ToInf_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
     	('dyInclusive50', '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
 	('Wjets', '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
      	('ttbar_lep', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('ttbar_lep50to500', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('ttbar_lep_500to800', '/TTToLL_MLL_500To800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('ttbar_lep_800to1200', '/TTToLL_MLL_800To1200_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('ttbar_lep_1200to1800', '/TTToLL_MLL_1200To1800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
	('ttbar_lep1800toInf', '/TTToLL_MLL_1800ToInf_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
	('Wantitop', '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
    	('tW', '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
# 
 	('qcd50to80', '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
    	('qcd80to120', '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd120to170', '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('qcd170to300', '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd300to470', '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('qcd470to600', '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd600to800', '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('qcd800to1000', '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd1000to1400', '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('qcd1400to1800', '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd1800to2400', '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
     	('qcd2400to3200', '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
 	('qcd3200', '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM'),

#	('RSGravitonM250','/RSGravToEEMuMu_kMpl-001_M-250_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM500','/RSGravToEEMuMu_kMpl-001_M-500_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM750','/RSGravToEEMuMu_kMpl-001_M-750_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM1000','/RSGravToEEMuMu_kMpl-001_M-1000_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM1500','/RSGravToEEMuMu_kMpl-001_M-1500_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM2000','/RSGravToEEMuMu_kMpl-001_M-2000_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM2500','/RSGravToEEMuMu_kMpl-001_M-2500_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM3000','/RSGravToEEMuMu_kMpl-001_M-3000_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM3500','/RSGravToEEMuMu_kMpl-001_M-3500_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#	('RSGravitonM4000','/RSGravToEEMuMu_kMpl-001_M-4000_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),

#	('CITo2Mu_Lam10TeVConLL', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_LLConM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
#	('CITo2Mu_Lam10TeVConLR', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_LRConM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
#	('CITo2Mu_Lam10TeVConRR', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_RRConM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
#	('CITo2Mu_Lam10TeVDesLL', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_LLDesM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
#	('CITo2Mu_Lam10TeVDesLR', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_LRDesM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
#	('CITo2Mu_Lam10TeVDesRR', '/CITo2Mu_GENSIM/szaleski-MuMu_miniAODSIM_RRDesM300_Use-28028af67189b3de7224b79195bd0e1d/USER'),
	
#	('DYTo2Mu_M300', '/DYTo2Mu_M300_CUETP8M1_13TeV_Pythia8_Corrected-v3/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),


	
        ]

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
    pobj = process.EventCounter * process.dileptonPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)


	
    if 'VBTF' not in cut_name and cut_name != 'Simple':
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
	for dataFilter in goodDataFiltersMiniAOD:
		#setattr(process,dataFilter 
		pobj = dataFilter * pobj


    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
        pobj = process.PrescaleToCommonMiniAOD * pobj ####### Now it seams that there are no prescaled path ########
    path = cms.Path(pobj)
    setattr(process, pathname, path)


def ntuplify(process, fill_gen_info=False):


    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    obj = process.prunedMCLeptons
    obj.src = cms.InputTag('prunedGenParticles')

    process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
                                           dimu_src = cms.InputTag('Our2016MuonsPlusMuonsMinus'),
					   met_src = cms.InputTag("slimmedMETs"),
					   jet_src = cms.InputTag("slimmedJets"),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					   #TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
 					   TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
                                           genEventInfo = cms.untracked.InputTag('generator'),
                                           metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter"))
                                           )
 
    if Electrons:
	    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction
        
    if hasattr(process, 'pathOur2016'):
	if fill_gen_info:
            process.pathOur2016 *=obj * process.SimpleNtupler 
	    if Electrons:
		process.pathOur2016 *=obj * process.SimpleNtupler * process.SimpleNtuplerEmu

ntuplify(process) #to have ntuples also running in interactive way

def for_mc(process, reco_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)
    switch_reco_process_name(process, reco_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

if isMC:
	for_mc(process,'PAT',True)


#def addGenFilter(process,name): 
#	print name

		

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

def for_data(process,GT):
    
    process.GlobalTag.globaltag = GT #RunH              #change line 52
    ntuplify(process)

if 'int_data' in sys.argv:
    for_data(process)
    printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
    printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

f = file('outfile_histos1', 'w')
f.write(process.dumpPython())
f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_datamc_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1
config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'global'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/jschulte'
config.Data.ignoreLocality = True 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_US_Purdue'
'''
    
    just_testing = 'testing' in sys.argv
        
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [

						('SingleMuonRun2016B-ReReco-v3', '/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
						('SingleMuonRun2016C-ReReco-v1', '/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
						('SingleMuonRun2016D-ReReco-v1', '/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
						('SingleMuonRun2016E-ReReco-v1', '/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
						('SingleMuonRun2016F-ReReco-v1', '/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
						('SingleMuonRun2016G-ReReco-v1', '/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD','80X_dataRun2_2016SeptRepro_v6'),
 						('SingleMuonRun2016H-PromptReco-v3', '/SingleMuon/Run2016H-PromptReco-v3/MINIAOD','80X_dataRun2_Prompt_v14'), #change global tag: 26 and 408
 						('SingleMuonRun2016H-PromptReco-v2', '/SingleMuon/Run2016H-PromptReco-v2/MINIAOD','80X_dataRun2_Prompt_v14'),  ##change global tag: 26 and 408

				# 							('SingleMuonRun2016G-PromptReco-v1',  '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD')

            ]

        lumi_lists = [
		# 'NoLumiMask'
		#           'DCSOnly',
		#            'Run2012PlusDCSOnlyMuonsOnly',
		'Run2016MuonsOnly',
		# 'Run2015',
		]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, GT, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name
	    print GT
            new_py = open('histosSimplified.py').read()
            new_py += "\nfor_data(process,'%s')\n"%GT
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
#config.Data.runRange = '256843-257490'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histosSimplified.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
           os.system('rm crabConfig.py histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 100000
    ''')



	ttbarFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
			       src = cms.InputTag('prunedGenParticles'),
			       min_mass = cms.double(50),
			       max_mass = cms.double(500), 
			       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	wwFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
				       src = cms.InputTag('prunedGenParticles'),
				       min_mass = cms.double(50),
				       max_mass = cms.double(200), 
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	     
	dyFilter = '''
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
				       src = cms.InputTag('prunedGenParticles'),                                      
				       )
for path_name, path in process.paths.iteritems():
	getattr(process,path_name).insert(0,process.DYGenMassFilter)'''

	



       
        for name, ana_dataset in samples:
            print name

            new_py = open('histosSimplified.py').read()
	    if "dyInclusive" in name:
		new_py += dyFilter
	    elif "ttbar_lep50to500" in name:
		new_py += ttbarFilter
	    elif "WWinclusive" in name:
		new_py += wwFilter
            new_py += "\nfor_mc(process,'PAT',True)\n"
            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histosSimplified.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

        if not just_testing:
		os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')

