#!/bin/tcsh
cmsRun -j ana_datamc_data.fjr.xml histos.py data >&! /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
gzip /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
mv /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out.gz .
fjr2json.py --output=ana_datamc_data.forlumi.json ana_datamc_data.fjr.xml
lumiCalc.py -i ana_datamc_data.forlumi.json overview > ana_datamc_data.lumi
mkdir ana_datamc_40ipb
mv ana_datamc_data.* ana_datamc_40ipb/

cmsRun -j ana_datamc_data.fjr.xml histos.py data all_good >&! /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
gzip /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out
mv /uscmst1b_scratch/lpc1/3DayLifetime/tucker/ana_datamc_data.out.gz .
fjr2json.py --output=ana_datamc_data.forlumi.json ana_datamc_data.fjr.xml
lumiCalc.py -i ana_datamc_data.forlumi.json overview > ana_datamc_data.lumi
mkdir ana_datamc_34ipb
mv ana_datamc_data.* ana_datamc_34ipb/
