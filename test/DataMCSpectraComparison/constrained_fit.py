from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import *
set_zp2mu_style()

ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)

ps = plot_saver('plots/constrained_fit', log=False, size=(600,600))

f = ROOT.TFile('data/Run2011MuonsOnly/ana_datamc_data.root')
t = f.SimpleNtupler.Get('t')

n = t.Draw('vertex_m:dil_mass:vertex_m_err', 'OurSelNew && dil_mass > 200')
g = ROOT.TGraph(n, t.GetV2(), t.GetV1())
g.SetMarkerStyle(7)
g.SetTitle(';unconstrained dimuon mass (GeV);dimuon vertex-constrained mass (GeV)')
g.Draw('AP')
g.GetXaxis().SetRangeUser(150, 1300)
g.GetYaxis().SetRangeUser(150, 1300)
g.GetYaxis().SetTitleOffset(1.2)
g.GetYaxis().SetLabelSize(0.03)
l = ROOT.TLine(150,150,1300,1300)
l.Draw()
ps.save('data_constrained_vs_un')

'''
ells = []
for i in xrange(n):
    m = dil_mass[i]
    v = vertex_m[i]
    e = vertex_m_err[i]
    ell = ROOT.TEllipse(m,v,e,e)
    ell.SetFillColor(17)
    ell.Draw()
    ells.append(ell)
g.Draw('P same')
l.Draw()
ps.save(name + '_constrained_vs_un_errors')
'''

f = ROOT.TFile('ana_datamc_Run2011AMuonsOnly/mc/ana_datamc_zssm1000.root')
t = f.SimpleNtupler.Get('t')

t.Draw('vertex_m/gen_dil_mass-1:dil_mass/gen_dil_mass-1', 'OurSelNew && abs(gen_dil_mass - 1000) < 20')
g = ROOT.TGraph(n, t.GetV2(), t.GetV1())
g.SetMarkerStyle(7)
g.Draw('AP')
g.GetXaxis().SetRangeUser(-0.5, 0.5)
g.GetYaxis().SetRangeUser(-0.5, 0.5)
g.GetYaxis().SetTitleOffset(1.2)
g.GetYaxis().SetLabelSize(0.03)
g.SetTitle(';unconstrained dimuon mass resolution (GeV);dimuon vertex-constrained mass resolution (GeV)')
ps.save('resolution_comparison')
