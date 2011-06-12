from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT, cumulative_histogram
from samples import *

ROOT.TH1.AddDirectory(False)

mc_fn_base = 'ana_datamc_current/muonsonly/mc/ana_datamc_%s.root'
rebin_factor = 5
int_lumi = 39.96

# mu+mu-, differential and cumulative
for cumulative in (False, True):
    histos = {}
    for sample in samples:
        f = ROOT.TFile(mc_fn_base % sample.name)
        h = f.OurMuonsPlusMuonsMinusHistos.Get('DileptonMass').Clone()
        h.Rebin(rebin_factor)
        h.Scale(sample.partial_weight * int_lumi)
        if cumulative:
            h = cumulative_histogram(h)
        histos[sample.name] = h

    f = ROOT.TFile('ana_datamc_current/muonsonly/ana_datamc_data.root')
    dataHist = f.OurMuonsPlusMuonsMinusHistos.Get('DileptonMass').Clone('dataHist')
    dataHist.Rebin(rebin_factor)
    if cumulative:
        dataHist = cumulative_histogram(dataHist)
        dataHist.SetName('dataHist')
    
    qcdHist = histos['inclmu15'].Clone('qcdHist')
    qcdHist.Add(histos['wjets'])

    ttbarHist = histos['ttbar'].Clone('ttbarHist')
    ttbarHist.Add(histos['singletop_tW'])
    ttbarHist.Add(histos['ww'])
    ttbarHist.Add(histos['wz'])
    ttbarHist.Add(histos['zz'])
    ttbarHist.Add(histos['ztautau'])
    
    zeeHist = histos['zmumu'].Clone('zeeHist')
    zeeHist.Add(histos['dy200'])
    zeeHist.Add(histos['dy500'])
    zeeHist.Add(histos['dy800'])

    zprime = histos['zssm1000'].Clone('zprime')
    
    # simulate stacking
    ttbarHist.Add(qcdHist)
    zeeHist.Add(qcdHist)
    zeeHist.Add(ttbarHist)

    fexport = ROOT.TFile('histos_export%s.root' % ('_cumulative' if cumulative else ''), 'RECREATE')
    dataHist.Write()
    zeeHist.Write()
    ttbarHist.Write()
    qcdHist.Write()
    zprime.Write()
    fexport.Close()

# e-mu
rebin_factor = 20
int_lumi = 35.47
histos = {}
for sample in samples:
    f = ROOT.TFile(mc_fn_base % sample.name)
    h = f.OurMuonsElectronsOppSignHistos.Get('DileptonMass').Clone()
    h.Rebin(rebin_factor)
    h.Scale(sample.partial_weight * int_lumi)
    histos[sample.name] = h

f = ROOT.TFile('ana_datamc_current/allgood/ana_datamc_data.root')
dataHist = f.OurMuonsElectronsOppSignHistos.Get('DileptonMass').Clone('dataHist')
dataHist.Rebin(rebin_factor)

jetsHist = histos['inclmu15'].Clone('jetsHist')
jetsHist.Add(histos['wjets'])

promptHist = histos['ttbar'].Clone('promptHist')

for x in ['zmumu', 'singletop_tW', 'ww', 'wz', 'zz', 'ztautau', 'dy200', 'dy500', 'dy800']:
    promptHist.Add(histos[x])
promptHist.Add(jetsHist)

fexport = ROOT.TFile('histos_export_emu.root', 'RECREATE')
dataHist.Write()
jetsHist.Write()
promptHist.Write()
fexport.Close()
