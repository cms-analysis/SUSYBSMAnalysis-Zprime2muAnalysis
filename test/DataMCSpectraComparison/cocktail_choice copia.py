from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ps = plot_saver('plots/cocktail_choice', log=False)

fdata = ROOT.TFile('ana_datamc_Run2011AMuonsOnly/ana_datamc_data.root')
fmc   = ROOT.TFile('ana_datamc_Run2011AMuonsOnly/mc/ana_datamc_dy200.root')

f = fdata

t = f.SimpleNtupler.Get('t')

ROOT.gStyle.SetOptStat(10)

def doit(name, cut):
    t.Draw('lep_cocktail_choice>>%s(5,0,5)' % name, 'OurSelNew && dil_mass > 400 && ' + cut)
    h = getattr(ROOT, name)
    h.Scale(1./h.GetEntries())
    h.SetTitle('Cocktail choice in %s;algo;fraction of muons passing our selection' % name)
    for i,n in enumerate(['global', 'tkonly', 'stalone', 'tpfms', 'picky']):
        h.GetXaxis().SetBinLabel(i+1, n)
    h.GetXaxis().LabelsOption('u')
    h.Draw('hist text00')
    ps.c.Update()
    s = h.GetListOfFunctions().FindObject('stats')
    s.SetX1NDC(0.14)
    s.SetY1NDC(0.78)
    s.SetX2NDC(0.40)
    s.SetY2NDC(0.85)
    ps.save(name)

doit('barrel',  'abs(lep_eta) < 0.85')
doit('overlap', '0.85 < abs(lep_eta) < 1.2')
doit('endcap',  '1.2 < abs(lep_eta)')
